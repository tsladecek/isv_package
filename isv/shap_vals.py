import shap
import numpy as np
import pandas as pd

from isv.config import settings
from isv.scripts.prepare_df import prepare
from isv.scripts.open_model import open_model
from isv.scripts.constants import HUMAN_READABLE, LOSS_ATTRIBUTES, GAIN_ATTRIBUTES


def shap_values_with_same_cnv_type(annotated_cnvs: pd.DataFrame, cnv_type: str):
    """Calculate SHAP values for CNVs with the same cnv type

    :param annotated_cnvs: Raw counts of genomic elements
    :param cnv_type: type of cnv
    :return: explainer object
    """
    X, X_train = prepare(annotated_cnvs, cnv_type, return_train=True)
    model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')

    explainer = shap.TreeExplainer(
        model,
        X_train,
        model_output='probability')

    return explainer(X).values


def shap_values(annotated_cnvs: pd.DataFrame):
    """Calculate SHAP values

    :param annotated_cnvs:
    :return: explainer object
    """
    del_ind = np.where(annotated_cnvs.cnv_type == "DEL")[0]
    dup_ind = np.where(annotated_cnvs.cnv_type == "DUP")[0]

    res = []

    if del_ind.shape[0] > 0:
        attributes = LOSS_ATTRIBUTES
        hr_attributes = ['SHAP_' + HUMAN_READABLE[i].replace(' ', '_') for i in attributes]
        sv = shap_values_with_same_cnv_type(annotated_cnvs.iloc[del_ind], "loss")
        res.append(pd.DataFrame(sv, columns=hr_attributes))

    if dup_ind.shape[0] > 0:
        attributes = GAIN_ATTRIBUTES
        hr_attributes = ['SHAP_' + HUMAN_READABLE[i].replace(' ', '_') for i in attributes]
        sv = shap_values_with_same_cnv_type(annotated_cnvs.iloc[dup_ind], "gain")
        res.append(pd.DataFrame(sv, columns=hr_attributes))

    res = pd.concat(res)
    res.index = np.concatenate([del_ind, dup_ind])
    res = res.sort_index()

    return res
