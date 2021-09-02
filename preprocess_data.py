#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import gzip
import os
import json
from isv.config import settings


def saveres(f, obj):
    """Save results to a gzipped json file
    :param f: filepath
    :param obj: object to save
    """
    with gzip.open(f, 'w') as g:    
        g.write(json.dumps(obj).encode('utf-8'))


def chrom_dict(arr, dummies=3):
    """Create per chromosome dictionary from input array
    :param arr: input array
    :return: chromosome dictionary
    """
    cd = {i:[] for i in range(1, 25)}
    
    for i in range(1, 25):
        temp = arr[np.where(arr[:, 0] == i)[0]].tolist()
        if len(temp) > 0:
            cd[i] = temp
        else:
            cd[i] = [[-1618] * dummies]  # Just a placeholder
        
    return cd


# %% gencode
print("Preprocessing Genes")
gencode_genes = pd.read_csv(os.path.join(settings.data_dir, "raw", "gencode_annotsv.tsv.gz"), sep='\t', compression='gzip')

# Only valid chromosomes
gencode_genes = gencode_genes.loc[[i in settings.valid_chromosomes for i in gencode_genes.chromosome]]

# rename types that contain "pseudogene" to just pseudogene
gencode_genes.type = ["pseudogene" if "pseudogene" in i else i for i in gencode_genes.type]

# only selected gene types. Rest will renamed to "unevaluated"
gene_types = ["protein_coding", "pseudogene", "lncRNA", "miRNA", "rRNA", "snRNA"]
gencode_genes.type = [i if i in gene_types else "unevaluated" for i in gencode_genes.type]

# Translate chromosome ids
gencode_genes.chromosome = [settings.chromosome_dict[i] for i in gencode_genes.chromosome]

# Translate gene types
gencode_genes.type = [settings.gene_type_dict[i] for i in gencode_genes.type]

gencode_genes.drop("gene_name", axis=1, inplace=True)
gencode_genes = gencode_genes.fillna(0)
gencode_genes = gencode_genes.values.astype(np.int64)

# chromdict
gencode_genes = chrom_dict(gencode_genes)

saveres(os.path.join(settings.data_dir, "preprocessed", "gencode_genes.json.gz"), gencode_genes)


# %% HI genes
print("Preprocessing HI genes")
hi_genes = pd.read_csv(os.path.join(settings.data_dir, "raw", "hi_genes.tsv.gz"), sep='\t', compression='gzip')

# Translate chromosome ids
hi_genes.chromosome = [settings.chromosome_dict[i] for i in hi_genes.chromosome]

# chromdict
hi_genes = chrom_dict(hi_genes.values.astype(np.int64))

saveres(os.path.join(settings.data_dir, "preprocessed", "hi_genes.json.gz"), hi_genes)

# %% HI/TS regions
print("Preprocessing HI/TS regions")
hits_regions = pd.read_csv(os.path.join(settings.data_dir, "raw", "hits_regions.tsv.gz"), sep='\t', compression='gzip')

# Translate chromosome ids
hits_regions.chromosome = [settings.chromosome_dict[i] for i in hits_regions.chromosome]

# chromdict
hits_regions = chrom_dict(hits_regions.values.astype(np.int64), dummies=5)

saveres(os.path.join(settings.data_dir, "preprocessed", "hits_regions.json.gz"), hits_regions)

# %% regulatory
print("Preprocessing Regulatory elements")
regulatory = pd.read_csv(os.path.join(settings.data_dir, "raw", "regulatory.tsv.gz"), sep='\t', compression='gzip')

# None type is "curated"
regulatory = regulatory.fillna("curated")

# Only valid chromosomes
regulatory = regulatory.loc[[i in settings.valid_chromosomes for i in regulatory.chromosome]]

# Translate chromosome ids
regulatory.chromosome = [settings.chromosome_dict[i] for i in regulatory.chromosome]

# Translate gene types
regulatory.type = [settings.regulatory_type_dict[i] for i in regulatory.type]

regulatory.drop("id", axis=1, inplace=True)

regulatory = regulatory.values.astype(np.int64)

# chromdict
regulatory = chrom_dict(regulatory)

saveres(os.path.join(settings.data_dir, "preprocessed", "regulatory.json.gz"), regulatory)
