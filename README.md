# ISV package

Python package for easy prediction of pathogenicity of annotated Copy Number Variants (CNVs)

The package contains two functions:

### 1. isv.predict(X_raw, cnv_type)
- returns an array of probabilities

### 2. shap_values(X_raw, cnv_type)
- calculates shap values for given CNVs

---
Both functions assume that input dataframes contain counts of following genomic elements:
```     
[
    'gencode_genes',
    'protein_coding',
    'pseudogenes',
    'mirna',
    'lncrna',
    'rrna',
    'snrna',
    'morbid_genes',
    'disease_associated_genes',
    'hi_genes',  # ONLY FOR LOSSES
    'regions_HI',  # ONLY FOR LOSSES
    'regions_TS',  # ONLY FOR GAINS
    'regulatory',
    'regulatory_enhancer',
    'regulatory_open_chromatin_region',
    'regulatory_promoter',
    'regulatory_promoter_flanking_region',
    'regulatory_ctcf_binding_site',
    'regulatory_tf_binding_site',
    'regulatory_curated'
]
```