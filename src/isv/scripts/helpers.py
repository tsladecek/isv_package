import pandas as pd
import numpy as np


def check_cnvs_obj(cnvs):
    if isinstance(cnvs, list) or isinstance(cnvs, np.ndarray):
        cnvs = pd.DataFrame(cnvs)
    assert isinstance(cnvs, pd.core.frame.DataFrame), "Input should be either list, np.ndarray or pd.DataFrame"
    assert cnvs.shape[1] == 4, "Input should have 4 columns: chromosome, start (GRCh38), end (GRCh38), cnv_type"
    assert sum(i not in ["DUP", "DEL"] for i in cnvs.iloc[:, 3]) == 0, "only 'DEL' and 'DUP' cnv_type values allowed"
    return cnvs
