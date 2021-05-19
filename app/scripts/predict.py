import xgboost as xgb
import pathlib

from app.config import settings
from app.scripts.prepare_df import prepare
from app.scripts.open_model import open_model


root_dir = pathlib.Path(__file__).parent.absolute()


def predict(cnv_type,
            data_path,
            proba=True):
    """Return model predictions for a selected dataframe
`
    :param cnv_type: cnv type, either "loss" or "gain"
    :param data_path: path to the dataframe to be predicted
    :param proba: return probabilities
    :returns: yhat: predicted values
    """
    model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')

    X = prepare(cnv_type=cnv_type, data_path=data_path)

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