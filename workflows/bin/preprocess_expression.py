#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse
import joblib

if __name__ == '__main__':
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--psam", required=True)
    parser.add_argument("--expr", required=True)
    parser.add_argument("--pca", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", required=False,default = 10)
    
    args = parser.parse_args()
    
    table = pd.read_csv(args.pca,sep = '\t', index_col = 0)
    table.columns = [f'g{x}' for x in table.columns]
    
    pca_model = PCA(n_components = args.n_components)
    
    expr_table = pd.read_csv(args.expr, index_col = 0, sep = '\t')
    expr_table = expr_table[expr_table.index.isin(table.index)]
    pca_df = pd.DataFrame(pca_model.fit_transform(expr_table), index = expr_table.index, columns = [f'ePC{i+1}' for i in range(args.n_components)])
    table = pca_df.join(table)
    
    scaler = StandardScaler()
    table = pd.DataFrame(scaler.fit_transform(table), index = table.index, columns = table.columns)
    
    table['SEX'] = pd.read_csv(args.psam, sep = '\t', index_col = 0)['SEX'] - 1
    table['intercept'] = 1

    genes = sorted(set(pd.read_csv(args.genes, header = None)[0])&set(expr_table.columns))
    expr_table = expr_table.loc[:, genes]

    X = table.loc[expr_table.index].values
    Y = np.log2(expr_table.values + 1)
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    expr_resid = pd.DataFrame((Y - (X @ beta)), index = expr_table.index, columns = expr_table.columns)

    expr_resid.to_csv(args.output, sep = '\t')

    coef = pd.DataFrame(beta, index = table.columns, columns = expr_table.columns).T

    joblib.dump({'beta':coef, 'pca_model':pca_model, 'scaler':scaler}, "regress_out_covariates.pkl")