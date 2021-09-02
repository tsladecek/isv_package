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

The package contains three functions:

### 1. isv.annotate(cnvs)
- annotates cnvs provided in a list, np.array or pandas DataFrame format represented in 4 columns: `chromosome`, `start (grch38)`, `end (grch38)` and `cnv_type`
- Returns an annotated dataframe which can be used as an input to following two functions

### 2. isv.predict(X_raw, cnv_type, proba)
- returns an array of isv predictions

### 3. isv.shap_values(X_raw, cnv_type)
- calculates shap values for given CNVs
