import json
import gzip
from sklearn_json import from_json, from_dict
import xgboost as xgb


def open_model(model_path):
    """Open and return a model from json file

    :param model_path: path to the model

    :return: model
    """
    if model_path.endswith('gz'):
        with gzip.open(model_path, 'r') as f:
            model = f.read()
            model = json.loads(model.decode('utf-8'))
            model = from_dict(model)
        return model

    else:
        with open(model_path, 'r') as f:
            a = f.readline()

        if a.startswith('{"learner"'):
            model = xgb.Booster()
            model.load_model(model_path)
            return model

        else:
            model = from_json(model_path)
            return model
