from .predict import predict
from .shap_vals import shap_values
from .annotate import annotate

import pandas as pd


def isv(cnvs, proba: bool = True, shap: bool = False):
    """Annotate and Predict pathogenicity of CNVs

    :param cnvs: a list, np.array or pandas dataframe with 4 columns representing chromosome (eg, chr3),
    cnv start (grch38), cnv end (grch38) and cnv_type (DUP or DEL)
    :param proba: whether probabilities should be returned
    :param shap: whether shap values should be calculated
    :return: ISV output as a pandas dataframe
    """
    final = cnvs.copy()
    annotated = annotate(cnvs)

    # 1. Add predictions
    final["ISV"] = predict(annotated, proba)

    # 2. Add SHAP values
    if shap:
        svs = shap_values(annotated)

        final = pd.concat([final, svs], axis=1)

    return final
