#!/usr/bin/env python

import pandas as pd
import numpy as np
import pickle
import argparse
from pgenlib import PgenReader

from scipy.stats import kruskal

def load_pfile(pgen, psam, pvar, snp_list, pvar_index):
    pgen_reader = PgenReader(bytes(pgen, 'utf8'))
    n_variants = len(snp_list)
    n_samples = pgen_reader.get_raw_sample_ct()

    genotype_matrix = np.empty((n_variants, n_samples), dtype=np.int8)

    for i,var in enumerate(snp_list):
        buf = np.empty(n_samples, dtype=np.int32)
        try:
            variant_idx = pvar_index[var]
            pgen_reader.read(variant_idx, buf)
            genotype_matrix[i, :] = buf
        except KeyError:
            genotype_matrix[i, :] = -1

    genotype_matrix = pd.DataFrame(genotype_matrix, index = snp_list, columns = pd.read_csv(psam, sep = '\t')['#IID']).T.replace(-1,np.nan)
    return genotype_matrix

if __name__ == '__main__':
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--pindex", required=True, help = ".pkl file indexing the given .pvar variant order") 
    parser.add_argument("--expr", required=True)
    parser.add_argument("--eqtl-folder", required=True)
    parser.add_argument("--tissue", required=False, default = 'whole_blood')
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()
    
    
    with open(args.pindex, 'rb') as file:
        pvar_index = pickle.load(file)
        
    eqtls = pd.read_csv(f'{args.eqtl_folder}/{args.tissue}.tsv', sep = '\t', index_col = ['ENSG','SNP'])['slope']
    expr = pd.read_csv(args.expr, sep = '\t', index_col = 0)
    
    res = []
    gene_list = list(set(eqtls.index.get_level_values(0).unique())&set(expr.columns))
    
    for ensg in gene_list:
        snp_list = list(eqtls.loc[ensg].index)

        genotype_matrix = load_pfile(args.pgen,args.psam,args.pvar,snp_list,pvar_index)
        genotype_matrix = genotype_matrix[genotype_matrix.index.isin(expr.index)]
        allele_freq = (genotype_matrix.sum(skipna=False)/(2*genotype_matrix.shape[0]))
        missingness = genotype_matrix.isna().sum()/genotype_matrix.shape[0]

        p = {}
        y = expr.loc[genotype_matrix.index, ensg]
        for snp in genotype_matrix.dropna(axis = 1).columns:
            x = genotype_matrix[snp]
            x_groups = {c:list(x[x==c].index) for c in x.unique()}
            y_groups = {c:y.loc[x_groups[c]] for c in x_groups}
            samples = y_groups.values()
            try:
                p[snp] = kruskal(*samples).pvalue
            except ValueError:
                pass

        p = pd.Series(p)

        ensg_res = pd.concat([pd.DataFrame({'maf':allele_freq,'missing':missingness,'kruskal_wallis_p':p})], keys = [ensg])
        res.append(ensg_res)
    res = pd.concat(res)
    
    if args.output.endswith('.csv'):
        sep = ','
    else:
        sep = '\t'
    
    res.to_csv(args.output, sep = sep)