import pathlib
import sys

import numpy as np
import pandas as pd

filepath_list = str(pathlib.Path(__file__).parent.absolute()).split('/')
ind = filepath_list.index('tests')
sys.path.insert(1, '/'.join(filepath_list[:ind]))

from isv.alternative import alternative_data, alternative_model


def test_model_retrain():
    np.random.seed(1618)

    for cnv_type in ['loss', 'gain']:
        cnv_type_ = ['DEL', 'DUP'][(cnv_type == 'gain') * 1]

        # load data
        data_dir = 'isv/data/'
        train = pd.read_csv(data_dir + f'train_{cnv_type}.tsv.gz', sep='\t', compression='gzip')
        validation = pd.read_csv(data_dir + f'validation_{cnv_type}.tsv.gz', sep='\t', compression='gzip')

        # extract chromosome, start and end coordinates
        train_cnvs, train_Y = train.loc[:, ["chr", "start_hg38", "end_hg38"]], train.clinsig.values
        val_cnvs, val_Y = validation.loc[:, ["chr", "start_hg38", "end_hg38"]], validation.clinsig.values

        # setup the 4th column - cnv_type
        train_cnvs["cnv_type"] = cnv_type_
        val_cnvs["cnv_type"] = cnv_type_

        # randomly generate extra columns
        train_extra = {'extra': np.random.rand(len(train_cnvs))}
        val_extra = {'extra': np.random.rand(len(val_cnvs))}

        # prepare data
        train_dmat = alternative_data(train_cnvs, train_Y, train_extra, cnv_type)
        val_dmat = alternative_data(val_cnvs, val_Y, val_extra, cnv_type)

        # Train Model
        model = alternative_model(train_dmat, val_dmat)

        # predict
        val_preds = model.predict(val_dmat)
        val_preds = (val_preds > 0.5) * 1

        assert len(val_preds) == len(val_Y)
