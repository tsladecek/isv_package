import xgboost as xgb
import shap

from isv.config import settings
from isv.scripts.prepare_df import prepare
from isv.scripts.open_model import open_model


def predict(X_raw, cnv_type, proba=True):
    """Return model predictions for a selected dataframe
`
    :param X_raw: Raw counts of genomic elements
    :param cnv_type: type of cnv
    :param proba: return probabilities
    :return: yhat: predicted values
    """
    model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')
    X = prepare(X_raw, cnv_type)

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