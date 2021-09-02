import xgboost as xgb

from isv.config import settings
from isv.scripts.prepare_df import prepare
from isv.scripts.open_model import open_model

import numpy as np
import pandas as pd


def predict_with_same_cnv_type(annotated_cnvs: pd.DataFrame, cnv_type: str, proba: bool = True):
    """Return model predictions for a selected dataframe
`
    :param annotated_cnvs: Raw counts of genomic elements
    :param cnv_type: type of cnv
    :param proba: return probabilities
    :return: yhat: predicted values
    """
    model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')
    X = prepare(annotated_cnvs, cnv_type)

    if proba:
        if isinstance(model, xgb.core.Booster):
            X_dmat = xgb.DMatrix(X)
            yhat = model.predict(X_dmat)
        else:
            yhat = model.predict_proba(X)[:, 1]

    else:
        if isinstance(model, xgb.core.Booster):
            X_dmat = xgb.DMatrix(X)
            yhat = model.predict(X_dmat)
            yhat = (yhat > 0.5) * 1
        else:
            yhat = model.predict(X)

    return yhat


def predict(annotated_cnvs: pd.DataFrame, proba: bool = True):
    """Predict bulk of CNVs with different cnv types

    :param annotated_cnvs: Annotated CNVs
    :param proba: whether probabilities should be calculated
    :return: predictions
    """

    del_ind = np.where(annotated_cnvs.cnv_type == "DEL")[0]
    dup_ind = np.where(annotated_cnvs.cnv_type == "DUP")[0]

    yh = np.empty(len(annotated_cnvs), dtype=[np.int8, np.float64][proba * 1])
    if len(del_ind) > 0:
        yh[del_ind] = predict_with_same_cnv_type(annotated_cnvs.iloc[del_ind], "loss", proba)

    if len(dup_ind) > 0:
        yh[dup_ind] = predict_with_same_cnv_type(annotated_cnvs.iloc[dup_ind], "gain", proba)

    return yh
