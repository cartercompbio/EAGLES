#!/usr/bin/env python

import pandas as pd
import numpy as np
from pgenlib import PgenReader
from sklearn.decomposition import PCA, FastICA as ICA
import argparse

def load_pgen_data(pgen_path, psam_path, pvar_path):
    """
    Load genotype data from PGEN format files into Python.
    
    Parameters:
    -----------
    pgen_path : str
        Path to the .pgen file (genotype data);
        if pgen_path is base.pgen, base.psam and base.pvar must exist
    
    Returns:
    --------
    genotype_matrix: pandas.DataFrame (sample x variant)
        values are 0-2 corresponding to count of ALT allele
    """
    
    # Load variant information (.pvar file)
    with open(pvar_path, 'r') as file:
        l = [x.strip().split('\t') for x in file if not x.startswith('##')]

    variants = pd.DataFrame(l[1:], columns = l[0])['ID'].values

    
    samples = pd.read_csv(psam_path, sep = '\t').rename({'#IID':'IID'},axis = 1)['IID'].values
            
    # Initialize the PGEN reader
    try:
        pgen_reader = PgenReader(bytes(pgen_path, 'utf8'))
    except Exception as e:
        raise Exception(f"Error opening PGEN file: {e}")    
    # Load genotype matrix
    n_variants = pgen_reader.get_variant_ct()
    n_samples = pgen_reader.get_raw_sample_ct()
        
    # Initialize genotype array
    genotype_matrix = np.empty((n_variants, n_samples), dtype=np.int8)
    
    # Read genotypes variant by variant
    buf = np.empty(n_samples, dtype=np.int32)
    for variant_idx in range(n_variants):
        pgen_reader.read(variant_idx, buf)
        genotype_matrix[variant_idx, :] = buf
        
    genotype_matrix = 2 - pd.DataFrame(genotype_matrix, index = variants, columns = samples).T
        
    return genotype_matrix

def pca_transform(g, thres = 0.999):
    pca = PCA()
    pc_df = pd.DataFrame(pca.fit_transform(g), index = g.index)
    
    pc_comp = pd.DataFrame(pca.components_, columns = g.columns, index = [f'PC{i}' for i in range(1,pca.components_.shape[0]+1)]).T

    pc_df = pc_df.loc[:, range((np.cumsum(pca.explained_variance_ratio_)<thres).sum())]
    pc_df.columns = [f'PC{i}' for i in range(1,pc_df.shape[1]+1)]
    return pc_df, pc_comp.loc[:, pc_df.columns]

def ica_transform(g, m_iter = 1000):
    ica = ICA(max_iter = m_iter)
    ic_df = pd.DataFrame(ica.fit_transform(g), index = g.index)

    ic_comp = pd.DataFrame(ica.components_, columns = g.columns, index = [f'IC{i}' for i in range(1,ica.components_.shape[0]+1)]).T

    ic_df.columns = [f'IC{i}' for i in range(1,ic_df.shape[1]+1)]
    return ic_df, ic_comp

def max_std_transform(g):
    g_df = g.loc[:, g.std().sort_values(ascending = False).index].copy()
    g_df.columns.names = [None]
    return g_df

def order_transform(g):
    order = sorted(g.columns, key = lambda x: int(x.split(':')[1]))
    return g.loc[:, order]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Variant feature engineering")
    parser.add_argument("--pgen", required=True, help="Path to pgen file")
    parser.add_argument("--psam", required=True, help="Path to psam file")
    parser.add_argument("--pvar", required=True, help="Path to pvar file")
    
    parser.add_argument("--method",
                        help="Transformation method: 'pca', 'ica',TBA...")
    parser.add_argument("--thres", required = False, default = 0.999, help = "if using pca method, proportion of explained varianced to be retained")
    
    parser.add_argument("--output", required=True, help="Path to output gene scores file")

    args = parser.parse_args()

    # Load data
    genotype = load_pgen_data(args.pgen, args.psam, args.pvar)



    # Then proceed with feature engineering method selected
    if args.method == "pca":
        result,loadings = pca_transform(genotype, args.thres)
    elif args.method == "ica":
        result,loadings = ica_transform(genotype)
    else:
        result = genotype
        loadings = pd.DataFrame([])

    result.to_csv(args.output, sep="\t")
    loadings.to_csv(args.output.replace('.tsv', '_loadings.tsv'), sep = "\t")
