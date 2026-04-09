import xgboost as xgb


def open_model(model_path):
    """Open and return a model from json file

    :param model_path: path to the model

    :return: model
    """
    model = xgb.Booster()
    model.load_model(model_path)
    return model
