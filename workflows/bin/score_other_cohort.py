#!/usr/bin/env python

import os
import joblib
import pandas as pd
import numpy as np
import argparse
from pgenlib import PgenReader
import shap

from fit_model import AlleleCount

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--model-folder", required=True, help="path to folder with trained eagles models (.pkl)")
    parser.add_argument("--covariates", required=False, default = None, help="Path to covariate file")
    parser.add_argument("--score-outfile", required=True, help = "Path to output file with all model predictions for this cohort")
    parser.add_argument("--missing-count-outfile", required=False, default = None, help = "Path to file describing number of missing features in this cohort")
    parser.add_argument("--replace-nan", required=False, type=int, default = None, help = "If given, constant value used to replace missing genotype info")
    args = parser.parse_args()

    if args.covariates is None:
        cov_table = pd.read_csv(args.covariates, sep = '\t', index_col = 0)
        cov_features = set(cov_table)
    else:
        cov_table = None
        cov_features = set()


    model_list = [x for x in os.listdir(args.model_folder) if x.endswith('.pkl')]
    model_snps = {}
    for f in model_list:
        tis = f[:f.find('ENSG')-1]
        ensg = f[len(tis)+1:-4]
        model_dict = joblib.load(f'{args.model_folder}/{f}')
        model_snps[(tis,ensg)] = list(model_dict['model'].feature_names_in_)

    all_snp_index = {x:None for y in model_snps for x in model_snps[y]}

    with open(args.pvar, 'r') as file:
        h = '##'
        while h.startswith('##'):
            h = next(file)

        h = h.split('\t')
        pos_index = h.index('ID')

        for i,line in enumerate(file):
            var_id = line.strip().split('\t')[pos_index]
            try:
                all_snp_index[var_id]
                all_snp_index[var_id] = i
            except:
                pass


    scores = {}
    missing_count = {}

    for (tis,ensg) in model_snps.keys():
        pgen_reader = PgenReader(bytes(args.pgen, 'utf8'))

        n_variants = len(model_snps[(tis,ensg)])
        n_samples = pgen_reader.get_raw_sample_ct()

        genotype_matrix = np.empty((n_variants, n_samples), dtype=np.int8)

        missing_snps = 0
        buf = np.empty(n_samples, dtype=np.int32)
        for i,var in enumerate(model_snps[(tis,ensg)]):
            variant_idx = all_snp_index[var]
            if variant_idx is not None:
                pgen_reader.read(variant_idx, buf)
                genotype_matrix[i, :] = buf
            else:
                genotype_matrix[i, :] = -1

        genotype_matrix = pd.DataFrame(genotype_matrix, index = model_snps[(tis,ensg)], columns = pd.read_csv(args.psam, sep = '\t')['#IID']).T.replace(-1, np.nan)
        if args.replace_nan is not None:
            genotype_matrix = genotype_matrix.fillna(args.replace_nan)

        model_dict = joblib.load(f'{args.model_folder}/{tis}_{ensg}.pkl')
        model = model_dict['model']

        X = genotype_matrix
        if args.covariates is not None:
            X = X.join(pd.read_csv(args.covariates, sep = '\t', index_col = 0))

        missing = set(model_dict['feature_names']) - set(X.columns)
        missing_count[(tis,ensg)] = {'missing':sum([1 for x in model_snps[(tis,ensg)] if all_snp_index[x] is None]) + len(missing),
                                     'total':len(model_dict['feature_names'])
                                    }
        missing_count[(tis,ensg)]['missing%'] = missing_count[(tis,ensg)]['missing']/missing_count[(tis,ensg)]['total']
        for col in missing:
            X[col] = np.nan

        X=X[model_dict['feature_names']]

        try:
            features = list(model_dict['scaler'].feature_names_in_)
            X_subset = X[features].copy()
            X_subset.columns = features

            X = pd.DataFrame(model_dict['scaler'].transform(X_subset), index = X_subset.index, columns = X_subset.columns)
            del X_subset, features
        except KeyError:
            pass

        scores[(tis,ensg)] = pd.Series(model.predict(X), index = X.index)

    
    missing_count = pd.DataFrame(missing_count).T.sort_values(by = 'missing%')
    scores = pd.DataFrame(scores).T.loc[missing_count.index]
                                                              
    
    
    scores.to_csv(args.score_outfile, sep = '\t')
    if args.missing_count_outfile is not None:
        missing_count.to_csv(args.missing_count_outfile, sep = '\t')


if __name__ == "__main__":
    main()