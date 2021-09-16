# ISV package

Python **pip** package for easy prediction of pathogenicity Copy Number Variants (CNVs)

---
## Install
#### Install with `pip install isv`
This will also automatically install all required additional packages. Thus it is recommended to install the package in a separate environment (e.g. virtualenv, conda, ...)

#### Package url: https://pypi.org/project/isv/

#### Module reference available at https://tsladecek.github.io/isv_package/

---
## Modules
##### The package contains a wrapper function:
### `isv.isv(cnvs, proba, shap)`
which automatically annotates and predicts `cnvs` provided in a list, np.array or pandas DataFrame format represented in 4 columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`

- The `proba` parameter controls whether probabilities should be calculated
- The `shap` parameter controls whether shap values should be calculated

#### and a Wrapper class (which is recommended):
### `isv.ISV(cnvs)`

with methods:
- ISV.predict(proba)
- ISV.shap(data=None)
  - where the `data` argument is optional
- ISV.waterfall(cnv_index)
  - for creating an interactive waterfall plot for a CNV at index `cnv_index`

---
#### The main subfunctions of the package are:

### 1. `isv.annotate(cnvs)`
- annotates cnvs provided in a list, np.array or pandas DataFrame format represented in 4 columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`
- Returns an annotated dataframe which can be used as an input to following two functions

### 2. `isv.predict(annotated_cnvs, proba)`
- returns an array of isv predictions. `annotated_cnvs` represents annotated cnvs returned by the annotate function

### 3. `isv.shap_values(annotated_cnvs)`
- calculates shap values for given CNVs. `annotated_cnvs` represents annotated cnvs returned by the annotate function

#### For example
1. using the simple wrapper
```
from isv import isv


cnvs = [
    ["chr8", 100000, 500000, "DEL"],
    ["chrX", 52000000, 55000000, "DUP"]
] 

results = isv(cnvs, proba=True, shap=True)
```

2. using the ISV class
```
from isv import ISV


cnvs = [
    ["chr8", 100000, 500000, "DEL"],
    ["chrX", 52000000, 55000000, "DUP"]
] 

cnv_isv = ISV(cnvs)
predictions = cnv_isv.predict(proba=True)
shap_vals = cnv_isv.shap()
cnv_isv.waterfall(cnv_index=1)
```

---
## Can be also used as a command line tool. Make sure to:

#### 1. clone the repository (https://github.com/tsladecek/isv_package)
#### 2. install requirements, e.g.

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### 3. Use ISV!
```
python isv_cmd.py -i <input_cnvs>.bed -o <outputpath> [-p] [-sv]
```
where the input should be a list of CNVs in a bed format, with columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`

Results will be saved in a tab separated file at path specified by user

Optionally, use following flags:
- **-p**: whether probabilities should be returned
- **-sv**: whether shap values should be calculated

#### For example

```
python isv_cmd.py -i examples/loss_gain_cnvs.bed -o examples/loss_gain_cnvs_out.bed -p -sv
```
