#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import numba as nb
import pandas as pd
import os
import gzip
import json
import time

from isv.config import settings


def open_data(f):
    """Open preprocessed data
    :param f: filepath
    :return: python dictionary of numpy arrays
    """
    with gzip.open(f, 'r') as g:
        temp = json.loads(g.read().decode("utf-8"))
    
    cd = {i:None for i in range(1, 25)} # just a placeholder 
    # as np arrays -> for numba
    for i in range(1, 25):
        cd[i] = np.array(temp[str(i)], dtype=np.int64)
    
    return cd


@nb.jit(nopython=True)
def get_el(arr, start: int, end: int):
    """Get overlapped elements
    :param arr: array to query
    :param start: start position
    :param end: end position
    :return: rows of array overlapped by start, end positions
    """
    return arr[np.where((arr[:, 1] <= end) & (arr[:, 2] >= start))[0]]


def annotate_cnv(chrom, start, end, gencode_genes, regulatory, hi_genes, hits_regions):
    """Annotate a candidate CNV

    :param chrom: chromosome identifier, eg. "chr1"
    :param start: start position on the GRCh38 assembly
    :param end: end position on the GRCh38 assembly
    :param gencode_genes: preprocessed gencode genes dict
    :param regulatory: preprocessed regulatory genes dict
    :param hi_genes: preprocessed hi_genes genes dict
    :param hi_regions: preprocessed hi_regions genes dict
    :param ts_regions: preprocessed ts_regions genes dict
    :return: annotation
    """
    
    overlapped_genes = get_el(gencode_genes[chrom], start, end)
    overlapped_regulatory = get_el(regulatory[chrom], start, end)
    overlapped_hi_genes = get_el(hi_genes[chrom], start, end)
    overlapped_hits_regions = get_el(hits_regions[chrom], start, end)
    
    return np.array([
        len(overlapped_genes),  # gencode_genes
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["protein_coding"]),
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["pseudogene"]),
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["miRNA"]),
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["lncRNA"]),
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["rRNA"]),
        np.sum(overlapped_genes[:, 3] == settings.gene_type_dict["snRNA"]),
        np.sum(overlapped_genes[:, 4]),  # morbid genes
        np.sum(overlapped_genes[:, 5]),  # disease associated genes
        len(overlapped_hi_genes),
        np.sum(overlapped_hits_regions[:, 3]),  # HI regions
        np.sum(overlapped_hits_regions[:, 4]),  # TS regions
        len(overlapped_regulatory),  # all overlapped regulatory elements
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["enhancer"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["open_chromatin_region"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["promoter"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["promoter_flanking_region"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["CTCF_binding_site"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["TF_binding_site"]),
        np.sum(overlapped_regulatory[:, 3] == settings.regulatory_type_dict["curated"])
    ], dtype=np.int32)


# %%
def annotate(cnvs):
    """
    :param cnvs: a list, np.array or pandas dataframe with 4 columns representing chromosome (eg, chr3),
    cnv start (grch38), cnv end (grch38) and cnv_type (DUP or DEL)
    :return: pd DataFrame of annotated CNVs
    """
    if isinstance(cnvs, list) or isinstance(cnvs, np.ndarray):
        cnvs = pd.DataFrame(cnvs)
    
    assert isinstance(cnvs, pd.core.frame.DataFrame),\
        "Please supply a list, np.array or pandas dataframe with 4 columns\
            representing the chromosome id, start (grch38), end (grch38)\
                and cnv_type (eiter DUP or DEL)"
    
    assert cnvs.shape[1] == 4, "Please supply dataframe with 4 columns only"

    cnvs.columns = ["chrom", "start", "end", "cnv_type"]

    assert np.sum([i in ["DEL", "DUP"] for i in cnvs.cnv_type]) == len(cnvs), "cnv type has to be either DUP or DEL"

    # Just in case something is wrong with indexes
    cnvs.reset_index(inplace=True, drop=True)

    # Make sure that chromosomes are in the right format
    cnvs.chrom = [i if i.startswith("chr") else f"chr{i}" for i in cnvs.chrom]

    # Load databases
    print("Loading databases")
    gencode_genes = open_data(os.path.join(settings.data_dir, "preprocessed", "gencode_genes.json.gz"))
    regulatory = open_data(os.path.join(settings.data_dir, "preprocessed", "regulatory.json.gz"))
    hi_genes = open_data(os.path.join(settings.data_dir, "preprocessed", "hi_genes.json.gz"))
    hits_regions = open_data(os.path.join(settings.data_dir, "preprocessed", "hits_regions.json.gz"))
    
    # annotate cnvs
    print("Annotating CNVs")
    annotated = np.empty((cnvs.shape[0], len(settings.attributes)), dtype=np.int64)
    start = time.time()
    for i, cnv in cnvs.iterrows():
        annotated[i] = annotate_cnv(settings.chromosome_dict[cnv.chrom],
                                    cnv.start, cnv.end, gencode_genes,
                                    regulatory, hi_genes, hits_regions)
    
    elapsed = round(time.time() - start, 9)
    print(f"Annotated in {elapsed} seconds")
    
    return pd.concat([cnvs, pd.DataFrame(annotated, columns=settings.attributes)], axis=1)
