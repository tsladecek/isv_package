import xgboost as xgb

from isv.config import settings
from isv.scripts.prepare_df import prepare
from isv.scripts.open_model import open_model

import numpy as np
import pandas as pd
import os


def predict_with_same_cnv_type(annotated_cnvs: pd.DataFrame, cnv_type: str):
    """Return model predictions for a selected dataframe
`
    :param annotated_cnvs: Raw counts of genomic elements
    :param cnv_type: type of cnv
    :return: yhat: predicted values
    """
    model = open_model(os.path.join(settings.model_dir, f'ISV_{cnv_type}.json'))
    X = prepare(annotated_cnvs, cnv_type)

    if isinstance(model, xgb.core.Booster):
        X_dmat = xgb.DMatrix(X)
        yhat = model.predict(X_dmat)
    else:
        yhat = model.predict_proba(X)[:, 1]

    return yhat


def predict(annotated_cnvs: pd.DataFrame, proba: bool = True, threshold: float = 0.95):
    """Predict bulk of CNVs with different cnv types

    :param annotated_cnvs: Annotated CNVs
    :param proba: whether probabilities should be calculated
    :param threshold: probability threshold for classifying CNVs into three classes: Pathogenic (>= threshold),
    Uncertain significance ((1-threshold, threshold)) or Benign (<= 1 - threshold)
    :return: predictions
    """

    del_ind = np.where(annotated_cnvs.cnv_type == "DEL")[0]
    dup_ind = np.where(annotated_cnvs.cnv_type == "DUP")[0]

    yh = np.empty(len(annotated_cnvs), dtype=np.float64)
    if len(del_ind) > 0:
        yh[del_ind] = predict_with_same_cnv_type(annotated_cnvs.iloc[del_ind], "loss")

    if len(dup_ind) > 0:
        yh[dup_ind] = predict_with_same_cnv_type(annotated_cnvs.iloc[dup_ind], "gain")

    if not proba:
        yh = np.array(["Pathogenic" if y >= threshold else "Benign" if y <= 1 - threshold \
                        else "Uncertain significance" for y in yh], dtype='O')

    return yh
