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


def chrom_dict(arr):
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
            cd[i] = [[1618, 1618, 1618]]  # Just a placeholder
        
    return cd


# %% annotsv genes
annotsv = pd.read_csv(os.path.join(settings.data_dir, "raw", "annotsv_genes.tsv"), sep='\t', low_memory=False)
annotsv = annotsv.iloc[np.where(annotsv.loc[:, "SV chrom"] != 'M')]

morbid = annotsv.query("morbidGenes == 'yes'").loc[:, "Gene name"].values
morbid = pd.DataFrame({"gene_name": morbid, "morbid": 1})

disease_assoc = annotsv.loc[:, "Gene name"][np.where(~annotsv.loc[:, "Mim Number"].isna())[0]].values
disease_assoc = pd.DataFrame({"gene_name": disease_assoc, "disease_associated": 1})

morbid_da = pd.merge(morbid, disease_assoc, "outer", "gene_name").drop_duplicates()

# %% gencode
gencode_gene_list = []

with gzip.open(os.path.join(settings.data_dir, "raw", "gencode.v37.annotation.gff3.gz"), 'r') as f:
    # Header
    for i in range(7):
        f.readline()
    # Data
    for line in f.readlines():
        line = line.decode('UTF-8')

        line_list = line.split('\t')

        if len(line_list) > 2 and line_list[2] == 'gene':
            chr_id, start, end = line_list[0], line_list[3], line_list[4]
            info = line_list[8].split(';')
            
            for i in info:
                if i.startswith('gene_name'):
                    gene_name = i.split('=')[1]

                if i.startswith('gene_type'):
                    gene_type = i.split('=')[1]

            gencode_gene_list.append([gene_name, chr_id, int(start), int(end), gene_type])

gencode_genes = pd.DataFrame(gencode_gene_list,
                             columns=['gene_name', 'chromosome', 'start', 'end', 'type'])


# Add morbid genes and disease associated genes
gencode_genes = pd.merge(gencode_genes, morbid_da, "left", "gene_name")

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

# %% ClinGen genes
clingen_genes = pd.read_csv(os.path.join(settings.data_dir, "raw", "ClinGen_gene_curation_list_GRCh38.tsv"),
                            sep='\t', skiprows=5)

hi_genes = clingen_genes.iloc[np.where(clingen_genes.loc[:, "Haploinsufficiency Score"] == 3)[0]]

hi_genes["chromosome"] = [i.split(':')[0] for i in hi_genes.loc[:, "Genomic Location"]]
hi_genes["start"] = [int(i.split(':')[1].split('-')[0].strip()) for i in hi_genes.loc[:, "Genomic Location"]]
hi_genes["end"] = [int(i.split(':')[1].split('-')[1].strip()) for i in hi_genes.loc[:, "Genomic Location"]]

hi_genes = hi_genes.loc[:, ["chromosome", "start", "end"]]

# Translate chromosome ids
hi_genes.chromosome = [settings.chromosome_dict[i] for i in hi_genes.chromosome]
hi_genes = hi_genes.values.astype(np.int64)

# chromdict
hi_genes = chrom_dict(hi_genes)

saveres(os.path.join(settings.data_dir, "preprocessed", "hi_genes.json.gz"), hi_genes)

# %% ClinGen HI regions
clingen_regions = pd.read_csv(os.path.join(settings.data_dir, "raw", "ClinGen_region_curation_list_GRCh38.tsv"),
                              sep='\t', skiprows=5)

hi_regions = clingen_regions.iloc[np.where(clingen_regions.loc[:, "Haploinsufficiency Score"] == 3)[0]]

hi_regions = hi_regions.loc[[i.startswith('chr') for i in hi_regions.loc[:, "Genomic Location"]]]

hi_regions["chromosome"] = [i.split(':')[0] for i in hi_regions.loc[:, "Genomic Location"]]
hi_regions["start"] = [int(i.split(':')[1].split('-')[0].strip()) for i in hi_regions.loc[:, "Genomic Location"]]
hi_regions["end"] = [int(i.split(':')[1].split('-')[1].strip()) for i in hi_regions.loc[:, "Genomic Location"]]

hi_regions = hi_regions.loc[:, ["chromosome", "start", "end"]]

# Translate chromosome ids
hi_regions.chromosome = [settings.chromosome_dict[i] for i in hi_regions.chromosome]
hi_regions = hi_regions.values.astype(np.int64)


# chromdict
hi_regions = chrom_dict(hi_regions)

saveres(os.path.join(settings.data_dir, "preprocessed", "hi_regions.json.gz"), hi_regions)

# %% ClinGen TS regions
clingen_regions = pd.read_csv(os.path.join(settings.data_dir, "raw", "ClinGen_region_curation_list_GRCh38.tsv"),
                              sep='\t', skiprows=5)

ts_regions = clingen_regions.iloc[np.where(clingen_regions.loc[:, "Triplosensitivity Score"] == 3)[0]]

ts_regions = ts_regions.loc[[i.startswith('chr') for i in ts_regions.loc[:, "Genomic Location"]]]

ts_regions["chromosome"] = [i.split(':')[0] for i in ts_regions.loc[:, "Genomic Location"]]
ts_regions["start"] = [int(i.split(':')[1].split('-')[0].strip()) for i in ts_regions.loc[:, "Genomic Location"]]
ts_regions["end"] = [int(i.split(':')[1].split('-')[1].strip()) for i in ts_regions.loc[:, "Genomic Location"]]

ts_regions = ts_regions.loc[:, ["chromosome", "start", "end"]]

# Translate chromosome ids
ts_regions.chromosome = [settings.chromosome_dict[i] for i in ts_regions.chromosome]
ts_regions = ts_regions.values.astype(np.int64)

# chromdict
ts_regions = chrom_dict(ts_regions)

saveres(os.path.join(settings.data_dir, "preprocessed", "ts_regions.json.gz"), ts_regions)

# %% regulatory
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
