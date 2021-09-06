from .predict import predict
from .shap_vals import shap_values
from .annotate import annotate
from .isv import ISV

import pandas as pd


def isv(cnvs, proba, shap):
    """Predict pathogenicity, and optionally calculate shap values of CNVs with this simple wrapper class

    :param cnvs: a list, np.array or pandas dataframe with 4 columns representing chromosome (eg, chr3),
    cnv start (grch38), cnv end (grch38) and cnv_type (DUP or DEL)
    :param proba: whether probabilities should be calculated
    :param shap: whether probabilities should be calculated
    :return: pandas dataframe of results
    """
    cnv_isv = ISV(cnvs)
    result = cnv_isv.predict(proba)
    if shap:
        temp = cnv_isv.shap()
        result = pd.concat([result, temp.iloc[:, 4:]], axis=1)

    return result
