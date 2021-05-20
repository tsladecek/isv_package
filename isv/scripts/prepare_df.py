import pandas as pd
from sklearn.preprocessing import RobustScaler

from isv.scripts.constants import LOSS_ATTRIBUTES, GAIN_ATTRIBUTES
from isv.config import settings


def prepare(X, cnv_type, return_train=False):
    """
    Extract relevant attributes for training and return training dataset
    together with labels, and scale the dataset - do same for validation dataset

    :param cnv_type: type of the cnv == ["loss", "gain"]
    :param X: pandas dataframe
    :param return_train: specify if transformed train data should be returned
    :return X: transformed dataframe if return_train is False. else tuple (X, X_train)
    """

    if cnv_type == 'loss':
        attributes = LOSS_ATTRIBUTES
    else:
        attributes = GAIN_ATTRIBUTES

    train_data_path = settings.data_dir + f'train_{cnv_type}.tsv.gz'
    X_train = pd.read_csv(train_data_path, compression='gzip', sep='\t')

    # Train
    X_train = X_train.loc[:, attributes]

    # evaluated data
    X_any = X.loc[:, attributes]

    # Scale
    scaler = RobustScaler()
    X_train = scaler.fit_transform(X_train)
    X_any = scaler.transform(X_any)

    if return_train:
        return X_any, X_train

    return X_any
