# ISV package

Python package for easy prediction of pathogenicity Copy Number Variants (CNVs)

---
## Requirements

```
xgboost>=1.4.2
shap>=0.39.0
sklearn-json>=0.1.0
numba>=0.53.1
numpy>=1.20.3
pandas>=1.2.4
numba>=0.53.1
```

The package might work with older versions, however above specified versions are recommended.
Make sure to install these packages before installing the `isv` package

---

The main function of the package is the function
### isv.isv(cnvs, proba, shap)
which automatically annotates and predicts `cnvs` provided in a list, np.array or pandas DataFrame format represented in 4 columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`

- The `proba` parameter controls whether probabilities should be calculated
- The `shap` parameter controls whether shap values should be calculated

And three main subfunction functions:

### 1. isv.annotate(cnvs)
- annotates cnvs provided in a list, np.array or pandas DataFrame format represented in 4 columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`
- Returns an annotated dataframe which can be used as an input to following two functions

### 2. isv.predict(annotated_cnvs, cnv_type, proba)
- returns an array of isv predictions. `annotated_cnvs` represents annotated cnvs returned by the annotate function

### 3. isv.shap_values(annotated_cnvs, cnv_type)
- calculates shap values for given CNVs. `annotated_cnvs` represents annotated cnvs returned by the annotate function

---
Can be also used as a command line tool. Make sure to:

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

File where results should be saved (as a tab-separated table)

Optionally, use following flags:
- **-p**: whether probabilities should be returned
- **-sv**: whether shap values should be calculated
