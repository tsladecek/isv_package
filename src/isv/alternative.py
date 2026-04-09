import sys
from typing import Union

import numpy as np
import pandas as pd
import xgboost as xgb

from isv import annotate
from isv.scripts.constants import GAIN_ATTRIBUTES, LOSS_ATTRIBUTES
from isv.scripts.helpers import check_cnvs_obj
from isv.scripts.prepare_df import prepare


def alternative_data(cnvs: Union[list, np.ndarray, pd.DataFrame], labels: Union[list, np.ndarray],
                     extra_columns: Union[pd.DataFrame, dict],
                     cnv_type: str) -> xgb.DMatrix:
    """Prepare dataset

    In the first step annotate and scale cnvs in bed format by ISV.
    Then add columns specified in a dictionary or dataframe.

    :param cnvs: cnvs specified in a bed format
    :param extra_columns: columns to be added to the dataset
    :param cnv_type: either DUP or DEL
    :return: pandas dataframe
    """
    cnvs = check_cnvs_obj(cnvs)

    if isinstance(extra_columns, dict):
        extra_columns = pd.DataFrame(extra_columns)

    cnv_type = cnv_type.lower()
    cnv_type_dict = {'del': 'DEL', 'loss': 'DEL', 'gain': 'DUP', 'dup': 'DUP'}
    try:
        cnv_type = cnv_type_dict[cnv_type]
    except KeyError:
        sys.exit("Only DEL/loss and DUP/gain values are allowed")

    cnv_type_set = set(cnvs.iloc[:, 3])
    assert len(cnv_type_set) == 1
    assert cnv_type in cnv_type_set

    annotated = annotate(cnvs)
    attributes = [LOSS_ATTRIBUTES, GAIN_ATTRIBUTES][cnv_type == 'DUP']
    annotated_scaled = pd.DataFrame(prepare(annotated, cnv_type), columns=attributes)

    # create dmat
    X = pd.concat([annotated_scaled, extra_columns], axis=1)
    y = np.array(labels).flatten()
    assert y.shape[0] == labels.shape[0]

    return xgb.DMatrix(X, y)


def alternative_model(train_dmat: xgb.DMatrix,
                      val_dmat: xgb.DMatrix,
                      params=None, num_boost_round: int = 100, early_stopping_rounds: int = 15, verbose_eval: int = 0):
    """Train a model with alternative

    :param train_dmat: result of "alternative_data" function for train cnvs
    :param val_dmat: result of "alternative_data" function for validation cnvs
    :param params: model parameters. If not set, ISV parameters are used
    :param num_boost_round: max number of boosting rounds
    :param early_stopping_rounds: Early stopping rounds
    :param verbose_eval: verbosity
    :return: xgboost model
    """
    train_Y = np.array(train_dmat.get_label())
    if params is None:
        params = {
            'max_depth': 8,
            'eta': 0.3,
            'gamma': 1,
            'subsample': 1,
            'lambda': 0.1,
            'colsample_bytree': 0.8,
            'scale_pos_weight': np.sqrt(sum(train_Y == 0) / sum(train_Y == 1)),
            'seed': 1618,
            'nthread': 4,
            'objective': 'binary:logistic',
            'eval_metric': 'logloss'}

    model = xgb.train(params, train_dmat, num_boost_round=num_boost_round, early_stopping_rounds=early_stopping_rounds,
                      evals=[(train_dmat, 'train'), (val_dmat, 'validation')], verbose_eval=verbose_eval)

    return model
