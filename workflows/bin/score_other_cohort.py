#!/usr/bin/env python

import os
import joblib
import pandas as pd
import numpy as np
import argparse
from pgenlib import PgenReader
import shap
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error as mse

from fit_model import AlleleCount

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--model-folder", required=True, help="path to folder with trained eagles models (.pkl)")
    parser.add_argument("--covariates", required=False, default = None, help="Path to covariate file")
    parser.add_argument("--outdir", required=True, help = "folder for output files")
    parser.add_argument("--replace-nan", required=False, type=int, default = None, help = "If given, constant value used to replace missing genotype info")
    parser.add_argument("--expr", required=False, help = "if given, will use to evaluate model predictions. Column names must match ENSG... of models in --model-folder")
    args = parser.parse_args()
    
    expr = pd.read_csv(args.expr, sep = '\t', index_col = 0)

    if args.covariates is not None:
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

    feature_summary = []
    model_performance = {}

    for (tis,ensg) in list(model_snps.keys()):
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
                missing_snps +=1

        genotype_matrix = pd.DataFrame(genotype_matrix, index = model_snps[(tis,ensg)], columns = pd.read_csv(args.psam, sep = '\t')['#IID']).T.replace(-1,np.nan)

        model_dict = joblib.load(f'{args.model_folder}/{tis}_{ensg}.pkl')
        model = model_dict['model']

        X = genotype_matrix
        if args.covariates is not None:
            X = X.join(pd.read_csv(args.covariates, sep = '\t', index_col = 0))

        missing = set(model_dict['feature_names']) - set(X.columns)

        for col in missing:
            X[col] = np.nan

        X=X[model_dict['feature_names']]

        if args.replace_nan is not None:
            X = X.fillna(args.replace_nan)

        try:
            features = list(model_dict['scaler'].feature_names_in_)
            X_subset = X[features].copy()
            X_subset.columns = features

            X = pd.DataFrame(model_dict['scaler'].transform(X_subset), index = X_subset.index, columns = X_subset.columns)
            del X_subset, features
        except KeyError:
            pass

        scores[(tis,ensg)] = pd.Series(model.predict(X), index = X.index)
        missing_count[(tis,ensg)] = {'snps':missing_snps,'total_snps':len(model_snps[(tis,ensg)]), 'missing_snps%':missing_snps/len(model_snps[(tis,ensg)])}

        # describe feature importances
        # for linear models this is the coefficient
        # otherwise this is absmax shap value within this cohort
        try:
            feat = pd.Series(model.coef_, index = model.feature_names_in_)
            feat = feat[(X.nunique()>1)&(feat!=0)]
            if feat.empty:
                best_snp = None
                snp_coef = None
            else:
                best_snp = feat.abs().idxmax()
                snp_coef = round(feat[best_snp]/abs(feat[best_snp]))
        except:
            explainer = shap.Explainer(model)
            shap_values = explainer(X)
            feat = pd.DataFrame(shap_values.values, index = X.index, columns = X.columns)
            feat = feat.abs().mean() #mean abs value used to prioritize important features for shap beeswarm plots...
            feat = feat[(X.nunique()>1)&(feat!=0)] #feature is relevant, and has variation within this cohort
            
            if feat.empty:
                best_snp = None
                snp_coef = None
            else:
                best_snp = feat.abs().idxmax()
                snp_coef = 2*(int(X[best_snp].corr(
                    pd.Series(shap_values.values[:, list(X.columns).index(best_snp)], index = X.index))>0) ) - 1

        if not feat.empty:
            feature_summary.append(pd.concat([pd.concat([feat], keys = [ensg])], keys = [tis]))

        if scores[(tis,ensg)].round(6).nunique()>1:
            try:
                y = expr[ensg]
                outlier = max([6*y.quantile(0.75) - 5*y.quantile(0.25),5])
                y = y[y<outlier]

                assert y.round(6).nunique()>1

                cur_scores = scores[(tis,ensg)][scores[(tis,ensg)].index.isin(y.index)]

                olap = list(set(X.index)&set(y.index)&set(cur_scores.index))
                cur_X = X.loc[olap]
                y = y.loc[olap]
                x2 = cur_scores.loc[olap]
                x1 = X.loc[olap, best_snp]*snp_coef

                x_df1 = pd.DataFrame({'x':x1})
                x_df2 = pd.DataFrame({'x':x2})

                m1 = LinearRegression().fit(x_df1,y)
                m2 = LinearRegression().fit(x_df2,y)

                model_performance[(tis,ensg)] = pd.Series({'r2_best_snp':x1.corr(y)**2,
                                                           'r2_score':x2.corr(y)**2, 
                                                           'r2_dif': x2.corr(y)**2 - x1.corr(y)**2,
                                                           'rmse_best_snp': mse(m1.predict(x_df1), y)**0.5,
                                                           'rmse_score': mse(m2.predict(x_df2), y)**0.5,
                                                           'rmse_score_minus_snp':mse(m2.predict(x_df2), y)**0.5 - mse(m1.predict(x_df1), y)**0.5
                                                          })
            except (AssertionError, TypeError):
                pass
    assert len(scores)>0
    if len(scores)>0:
        scores = pd.DataFrame(scores)
        scores.to_csv(f'{args.outdir}/scores.tsv', sep = '\t')
        
    if len(missing_count)>0:
        missing_count = pd.DataFrame(missing_count).T
        missing_count.to_csv(f'{args.outdir}/missing_counts.tsv', sep = '\t')
        
    if len(feature_summary)>0:    
        feature_summary = pd.concat(feature_summary).reset_index()
        feature_summary.columns = ['TISSUE','GENE','SNP','SHAP']
        feature_summary.to_csv(f'{args.outdir}/feature_summary.tsv', sep = '\t', index = False)
    
    if len(model_performance)>0:
        model_performance = pd.DataFrame(model_performance).T
        model_performance.to_csv(f'{args.outdir}/model_performance.tsv', sep = '\t')


if __name__ == "__main__":
    main()