#!/usr/bin/env python

import os
import joblib
import pandas as pd
import numpy as np
import argparse
from pgenlib import PgenReader
import shap
import statsmodels.api as sm
import pickle

from fit_model import AlleleCount
from sklearn.linear_model import LinearRegression

def index_cohort(pvar_in, pkl_out):
    snp_index = {}
    with open(pvar_in, 'r') as file:
        h = '##'
        while h.startswith('##'):
            h = next(file)

        h = h.split('\t')
        pos_index = h.index('ID')

        for i,line in enumerate(file):
            var_id = line.strip().split('\t')[pos_index]
            snp_index[var_id] = i

    with open(pkl_out, 'wb') as file:
        pickle.dump(snp_index, file)

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--pindex", required=True, help = ".pkl file indexing the given .pvar variant order") 
    parser.add_argument("--model-folder", required=True, help="path to folder with trained eagles models (.pkl)")
    parser.add_argument("--outdir", required=True, help = "folder for output files")
    
    parser.add_argument("--covariates", required=False, default = None, help="Path to covariate file. If not given, numeric columns in psam file will be used for covariates")
    parser.add_argument("--replace-nan", required=False, type=int, default = None, help = "If given, constant value used to replace missing genotype info")
    parser.add_argument("--expr", required=False, default=None,help = "if given, will use to evaluate model predictions. Column names must match ENSG... of models in --model-folder")
    args = parser.parse_args()
    
    if args.expr is not None:
        expr = pd.read_csv(args.expr, sep = '\t', index_col = 0)
    else:
        expr = pd.DataFrame([])

    if args.covariates is not None:
        cov_table = pd.read_csv(args.covariates, sep = '\t', index_col = 0)
    else:
        cov_table = pd.read_csv(args.psam, sep = '\t', index_col = 0)
        cov_table = cov_table.select_dtypes(include='number')
        
    cov_features = set(cov_table)
    
    model_paths = {(f[:f.find('ENSG')-1], f[f.rfind('_')+1:-4]):f for f in os.listdir(args.model_folder) if f.endswith('.pkl')}
    with open(args.pindex, 'rb') as file:
        pvar_index = pickle.load(file)
        

    scores = {}
    missing_count = {}
    feature_summary = []
    model_performance = {}

    for (tis,ensg) in model_paths.keys():
        model_dict = joblib.load(f'{args.model_folder}/{tis}_{ensg}.pkl')
        model = model_dict['model']
        
        
        ###
        # attempt to load all features from plink files
        # count missing features
        ###
        
        pgen_reader = PgenReader(bytes(args.pgen, 'utf8'))
        
        n_variants = len(model.feature_names_in_)
        n_samples = pgen_reader.get_raw_sample_ct()

        genotype_matrix = np.empty((n_variants, n_samples), dtype=np.int8)

        missing_ct = 0
        for i,var in enumerate(model.feature_names_in_):
            buf = np.empty(n_samples, dtype=np.int32)
            try:
                variant_idx = pvar_index[var]
                pgen_reader.read(variant_idx, buf)
                genotype_matrix[i, :] = buf
            except KeyError:
                genotype_matrix[i, :] = -1
                missing_ct +=1

        genotype_matrix = pd.DataFrame(genotype_matrix, index = model.feature_names_in_, columns = pd.read_csv(args.psam, sep = '\t')['#IID']).T.replace(-1,np.nan)

        ###
        # attempt to find missing features in the covariate table
        # if any are found, update missing count
        ###
        X = genotype_matrix
        for col in cov_features&set(model.feature_names_in_):
            X[col] = cov_table[col]
            missing_ct-=1
        X=X[model.feature_names_in_]
        
        if args.replace_nan is not None:
            X = X.fillna(args.replace_nan)
        elif isinstance(model, LinearRegression):
            X = X.fillna(0)
            
        ###
        # if features scaled during training, apply same scaler to X
        ###

        try:
            features = list(model_dict['scaler'].feature_names_in_)
            X_subset = X[features].copy()
            X_subset.columns = features

            X = pd.DataFrame(model_dict['scaler'].transform(X_subset), index = X_subset.index, columns = X_subset.columns)
            del X_subset, features
        except KeyError:
            pass

        scores[(tis,ensg)] = pd.Series(model.predict(X), index = X.index)
        missing_count[(tis,ensg)] = {'missing_features':missing_ct,'total_features':len(model.feature_names_in_), 'missing_feature%':missing_ct/len(model.feature_names_in_)}

        ###
        # describe feature importances
        # absmax shap value within this cohort
        # save to feature_summary
        ###
        try:
            explainer = shap.Explainer(model, X)
            shap_values = explainer(X)
        except:
            explainer = shap.TreeExplainer(model, feature_perturbation='tree_path_dependent' )
            shap_values = explainer(X)  

        feat = pd.DataFrame(shap_values.values, index = X.index, columns = X.columns)
        feat = feat.abs().mean()
        feat = feat[(X.nunique()>1)&(feat!=0)] #feature is relevant, and has variation within this cohort

        if feat.empty:
            best_snp = None
        else:
            best_snp = feat.idxmax()

        if not feat.empty:
            feature_summary.append(pd.concat([pd.concat([feat], keys = [ensg])], keys = [tis]))
            
        ###
        # if variation exists within current cohort's predictions
        # and variation exists within current cohort's phenotype
        # 
        # fit 4 linear models: y~x+1
        #     y: gene expression
        #     x: (1) allele count of best_snp (2) EAGLES predictions (3) best_snp,covariates (4) EAGLES predictions, covariates
        #
        # reports for each linear model
        #     r2
        #     rmse
        #     coef of EAGLES prediction or best_snp
        #     pvalue of EAGLES prediction or best_snp
        #     intercept


        if scores[(tis,ensg)].round(6).nunique()>1:
            try:
                y = expr[ensg]
                outlier = max([6*y.quantile(0.75) - 5*y.quantile(0.25),5])
                y = y[y<outlier]

                assert y.round(6).nunique()>1

                cur_scores = scores[(tis,ensg)]

                olap = list(set(X.index)&set(y.index)&set(cur_scores.index))

                y = y.loc[olap]
                x2 = cur_scores.loc[olap]
                x1 = X[best_snp].loc[olap]

                X1 = sm.add_constant(x1)
                model1 = sm.OLS(y, X1)
                results1 = model1.fit()

                X2 = sm.add_constant(x2)
                model2 = sm.OLS(y,X2)
                results2 = model2.fit()

                X3 = sm.add_constant(x1).join(cov_table)
                model3 = sm.OLS(y,X3)
                results3 = model3.fit()

                X4 = sm.add_constant(x2).join(cov_table)
                model4 = sm.OLS(y,X4)
                results4 = model4.fit()

                model_performance[(tis,ensg)] = pd.Series({'r2|EAGLES':results2.rsquared, 
                                                           'rmse|EAGLES': results2.mse_resid**0.5,
                                                           'coef|EAGLES': results2.params[0],
                                                           'intercept|EAGLES': results2.params['const'],
                                                           'pval|EAGLES': results2.pvalues[0],
                                                           
                                                           'r2|best_snp':results1.rsquared,
                                                           'rmse|best_snp': results1.mse_resid**0.5,
                                                           'coef|best_snp': results1.params[best_snp],
                                                           'intercept|best_snp':results1.params['const'],
                                                           'pval|best_snp':results1.pvalues[best_snp],
                                                           
                                                           'r2|EAGLES_plus_cov':results4.rsquared, 
                                                           'rmse|EAGLES_plus_cov': results4.mse_resid**0.5,
                                                           'coef|EAGLES_plus_cov': results4.params[0],
                                                           'intercept|EAGLES_plus_cov': results4.params['const'],
                                                           'pval|EAGLES_plus_cov': results4.pvalues[0],
                                                           
                                                           'r2|best_snp_plus_cov':results3.rsquared,
                                                           'rmse|best_snp_plus_cov': results3.mse_resid**0.5,
                                                           'coef|best_snp_plus_cov': results3.params[best_snp],
                                                           'intercept|best_snp_plus_cov':results3.params['const'],
                                                           'pval|best_snp_plus_cov':results3.pvalues[best_snp]
                                                        })
            except (AssertionError, TypeError, KeyError):
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
