import shap

from isv.config import settings
from isv.scripts.prepare_df import prepare
from isv.scripts.open_model import open_model


def shap_values(X_raw, cnv_type):
    """Calculate SHAP values

    :param X_raw: Raw counts of genomic elements
    :param cnv_type: type of cnv
    :return: explainer object
    """
    X, X_train = prepare(X_raw, cnv_type, return_train=True)
    model = open_model(settings.model_dir + f'ISV_{cnv_type}.json')

    explainer = shap.TreeExplainer(
        model,
        X_train,
        model_output='probability')

    return explainer(X)
