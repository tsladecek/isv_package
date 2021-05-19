import pandas as pd
from sklearn.preprocessing import RobustScaler

from app.scripts.constants import LOSS_ATTRIBUTES, GAIN_ATTRIBUTES
from app.config import settings


def prepare(cnv_type,
            data_path):
    """
    Extract relevant attributes for training and return training dataset
    together with labels, and scale the dataset - do same for validation dataset

    :param cnv_type: type of the cnv == ["loss", "gain"]
    :param data_path: path to the data to predict
    :return X: transformed dataframe
    """

    if cnv_type == 'loss':
        attributes = LOSS_ATTRIBUTES
    else:
        attributes = GAIN_ATTRIBUTES

    train_data_path = settings.data_dir + f'train_{cnv_type}.tsv.gz'
    X_train = pd.read_csv(train_data_path, compression='gzip', sep='\t')
    X_any = pd.read_csv(data_path, compression='gzip', sep='\t')

    # Train
    X_train = X_train.loc[:, attributes]

    # evaluated data data
    X_any = X_any.loc[:, attributes]

    # Scale
    scaler = RobustScaler()

    scaler.fit(X_train)

    X_any = scaler.transform(X_any)

    return X_any
