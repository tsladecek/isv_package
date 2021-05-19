import xgboost as xgb
import shap

from app.config import settings
from app.scripts.prepare_df import prepare
from app.scripts.open_model import open_model


class Predict:
    def __init__(self, cnv_type, data_path):
        """
        Predict class

        :param cnv_type: cnv type, either "loss" or "gain"
        :param data_path: path to the dataframe to be predicted
        """
        self.cnv_type = cnv_type
        self.model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')
        self.X = prepare(cnv_type=cnv_type, data_path=data_path)

    def predict(self,
                proba=True):
        """Return model predictions for a selected dataframe
    `
        :param proba: return probabilities
        :returns: yhat: predicted values
        """
        model = self.model
        X = self.X

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

    def shap_values(self):

        X_train = prepare(self.cnv_type, settings.data_dir + f'train_{self.cnv_type}.tsv.gz')

        explainer = shap.TreeExplainer(
            self.model,
            X_train,
            model_output='probability')
        return self.explainer(self.X)
