import sys
import pathlib
import pandas as pd

filepath_list = str(pathlib.Path(__file__).parent.absolute()).split('/')
ind = filepath_list.index('tests')
sys.path.insert(1, '/'.join(filepath_list[:ind]))

from isv import isv, ISV

cnvs = [['chrX', 50000, 10000, "DEL"], ["chr7", 50, 600000, "DUP"]]


def test_isv_wrapper_fun():
    res = isv(cnvs, proba=True, shap=True)

    assert isinstance(res, pd.core.frame.DataFrame)


def test_isv_class():
    isv_cnvs = ISV(cnvs)

    p = isv_cnvs.predict()
    s = isv_cnvs.shap()
