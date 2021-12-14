import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

from isv.annotate import annotate
from isv.predict import predict as predict_cnvs
from isv.scripts.constants import LOSS_ATTRIBUTES, GAIN_ATTRIBUTES, HUMAN_READABLE, DESCRIPTIONS
from isv.scripts.helpers import check_cnvs_obj
from isv.shap_vals import shap_values, shap_values_with_same_cnv_type


class ISV:
    """Annotate and Predict pathogenicity of CNVs

    :param cnvs: a list, np.array or pandas dataframe with 4 columns representing chromosome (eg, chr3), \
    cnv start (grch38), cnv end (grch38) and cnv_type (DUP or DEL)

    :return: ISV output as a pandas dataframe
    """

    def __init__(self, cnvs):

        cnvs = check_cnvs_obj(cnvs)

        cnvs.columns = ["chromosome", "start", "end", "cnv_type"]

        self.annotated = annotate(cnvs)
        self.cnvs = cnvs

    def predict(self, proba: bool = True, threshold: float = 0.95):
        """Generate ISV predictions

        :param proba: whether probabilities should be calculated
        :param threshold: probability threshold for classifying CNVs into three classes: Pathogenic (>= threshold), \
        Uncertain significance ((1-threshold, threshold)) or Benign (<= 1 - threshold)

        :return: dataframe with last column representing the ISV predictions
        """
        res = self.cnvs.copy()
        res["ISV"] = predict_cnvs(self.annotated, proba=proba, threshold=threshold)
        return res

    def shap(self, df: pd.core.frame.DataFrame = None):
        """Calculate SHAP values

        :return: dataframe of shap values
        """
        if df is not None:
            return shap_values(df)

        sv = shap_values(self.annotated)
        return pd.concat([self.cnvs, sv], axis=1)

    def waterfall(self,
                  cnv_index: int,
                  filepath: str = "temp-plot.html",
                  return_fig: bool = False,
                  pathogenic_color: str = "rgb(255, 0, 50)",
                  benign_color: str = "rgb(58, 130, 255)",
                  text_position: str = 'outside',  # 'none' for no text
                  width: int = 800,
                  height: int = 800,
                  ):
        """Waterfall plot for CNV at specified index

        :param cnv_index: CNV which should be plotted
        :param filepath: Path where resulting html file will be saved
        :param return_fig: whether raw figure should be returned
        :param pathogenic_color: color of bars pushing predictions to pathogenic values
        :param benign_color: color of bars pushing predictions to benign values
        :param text_position: text position
        :param width: figure width
        :param height: figure height

        :return: html plot
        """
        cnv_type = ["loss", "gain"][(self.cnvs.iloc[cnv_index:(cnv_index + 1)].cnv_type.item() == "gain") * 1]

        sv = shap_values_with_same_cnv_type(self.annotated.iloc[cnv_index:(cnv_index + 1)],
                                            cnv_type,
                                            raw=True)

        visibility_dict = {
            'default': [True, False, False],
            'sorted': [False, True, False],
            'abs_sorted': [False, False, True]
        }

        attributes = [LOSS_ATTRIBUTES, GAIN_ATTRIBUTES][(cnv_type == 'gain') * 1]
        feature_names = [HUMAN_READABLE[i] for i in attributes]
        descriptions = [DESCRIPTIONS[i] for i in attributes]

        # DataFrame for Bar Charts
        data = pd.DataFrame({'Feature': feature_names,  # Nice Feature without SHAP prefix
                             'Descriptions': descriptions,
                             'Value': sv.data[0],
                             'SHAP': sv.values[0],
                             'raw': self.annotated.loc[:, attributes].iloc[cnv_index].values,
                             'colors': [pathogenic_color if i > 0 else benign_color for i in sv.values[0]]})

        data = data.iloc[::-1]

        # sort data by SHAP value
        sorted_data = data.iloc[np.argsort(data.SHAP)]

        # sort data by absoulte SHAP value
        abs_sorted_data = data.iloc[np.argsort(np.abs(data.SHAP))]

        # hover text formating
        hover_text = '{}<br>Value: {:2.3f}'

        v = visibility_dict['sorted']  # default visibility

        fig = go.Figure()

        for i, temp in enumerate([data, sorted_data, abs_sorted_data]):
            # add value to the left of the name of the feature
            temp["alt_Feature"] = [
                f'<span style="font-size: 10px; color: gray">({temp.raw.iloc[i]})</span> = <span style="font-size: ' \
                f'14px;">{temp.Feature.iloc[i]}</span> '
                for i in range(len(temp))]

            # Main Bar Plot
            fig.add_trace(
                go.Waterfall(
                    x=temp.SHAP,
                    y=temp.alt_Feature,
                    base=sv.base_values[0],
                    orientation='h',
                    hovertext=[hover_text.format(temp.Descriptions.iloc[i], temp.raw.iloc[i]) for i in
                               range(len(temp))],
                    hoverinfo="text + delta + initial",
                    connector={"mode": "between", "line": {"width": 0.2, "color": "gray", "dash": "solid"}},
                    decreasing={
                        "marker": {"color": benign_color, "line": {'width': 0.1}}},
                    increasing={"marker": {"color": pathogenic_color, "line": {'width': 0}}},
                    name='Pathogenic/Benign',
                    visible=v[i],
                    text=[np.round(i, 2) if i <= 0 else f'+{np.round(i, 2)}' for i in temp.SHAP],
                    textposition=text_position,
                    textfont={'color': 'gray'}
                )
            )

        # Buttons
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    active=1,
                    buttons=list([
                        dict(
                            label="Default",
                            args=[
                                {"visible": visibility_dict['default']}
                            ]
                        ),
                        dict(
                            label="Sorted",
                            args=[
                                {"visible": visibility_dict['sorted']}
                            ],
                        ),
                        dict(
                            label="Abs Sorted",
                            args=[
                                {"visible": visibility_dict['abs_sorted']}
                            ]
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.11,
                    xanchor="left",
                    y=1.1,
                    yanchor="top"
                ),
            ]
        )

        fig.add_shape(type='line',
                      x0=sv.base_values[0],
                      y0=0,
                      x1=sv.base_values[0],
                      y1=len(data),
                      line=dict(color='grey', width=0.3, dash='dash'),
                      )

        # General Layout
        chrom, start, end, cnv_type = self.cnvs.iloc[cnv_index].values
        fig.update_layout(
            template='plotly_white',
            # showlegend=True,
            xaxis_title=None,
            yaxis_title=None,
            title=dict(text=f"{['Deletion', 'Duplication'][cnv_type == 'DUP']} of {chrom}:{start}-{end}"),
            legend=dict(
                itemclick=False,
                itemdoubleclick=False,
                orientation="h",
                yanchor="bottom",
                y=1,
                xanchor="right",
                x=0.5
            ),
            width=width,
            height=height,
            xaxis={"showgrid": True,
                   "nticks": 5,
                   "range": [sv.base_values[0] + np.min(np.cumsum(sv.values[0])) - 0.25, \
                             sv.base_values[0] + np.max(np.cumsum(sv.values[0])) + 0.1]},
            hoverlabel=dict(
                font_size=16,
                font_family="Rockwell",
                font=dict(color='white')
            ),
            margin=dict(t=20, b=20, l=180, r=0),
        )

        if return_fig:
            return fig
        plot(fig, filename=filepath)
