#!/usr/bin/env python

import os
import joblib
import pandas as pd
import numpy as np
import argparse
from pgenlib import PgenReader
import shap
import pickle

from sklearn.linear_model import LinearRegression
from sklearn.utils import get_tags
from sklearn.metrics import r2_score, root_mean_squared_error as rmse

import xgboost as xgb

import statsmodels.api as sm

def supports_nan(model):
    try:
        return get_tags(model).input_tags.allow_nan
    except Exception:
        return False 

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--pindex", required=True, help = ".pkl file indexing the given .pvar variant order") 
    parser.add_argument("--model-folder", required=True, help="path to folder with trained eagles models (.pkl)")
    parser.add_argument("--outdir", required=True, help = "folder for output files")    
    parser.add_argument("--expr", required=False, default=None,help = "if given, will use to evaluate model predictions. Column names must match ENSG... of models in --model-folder")
    args = parser.parse_args()
    
    if args.expr is not None:
        expr = pd.read_csv(args.expr, sep = '\t', index_col = 0)
    else:
        expr = pd.DataFrame([])
    
    model_paths = {(f[:f.find('ENSG')-1], f[f.rfind('_')+1:-4]):f for f in os.listdir(args.model_folder) if f.endswith('.pkl')}
    with open(args.pindex, 'rb') as file:
        pvar_index = pickle.load(file)

    try:
        columns = pd.read_csv(args.psam, sep = '\t')['#IID']
    except KeyError:
        columns = pd.read_csv(args.psam, sep = '\t')['IID']
    
    scores = {}
    missing_count = {}
    feature_summary = []
    model_performance = {}
    
    model_paths = {(f[:f.find('ENSG')-1], f[f.find('ENSG'):-4]):f for f in os.listdir(args.model_folder) if f.endswith('.pkl')}
    for (tis,ensg) in model_paths.keys():
        
        model_dict = joblib.load(f'{args.model_folder}/{tis}_{ensg}.pkl')
        model = model_dict['model']
        
        pgen_reader = PgenReader(bytes(args.pgen, 'utf8'))
        
        n_variants = len(model.feature_names_in_)
        n_samples = pgen_reader.get_raw_sample_ct()
        
        X = np.empty((n_variants, n_samples), dtype=np.int8)
        
        missing_ct = 0
        for i,var in enumerate(model.feature_names_in_):
            buf = np.empty(n_samples, dtype=np.int32)
            try:
                variant_idx = pvar_index[var]
                pgen_reader.read(variant_idx, buf)
                X[i, :] = buf
            except KeyError:
                X[i, :] = -1
                missing_ct +=1
        

        X = pd.DataFrame(X, index = model.feature_names_in_, columns=columns).T.replace(-1,np.nan)
        geno = X
        
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
        
        ###
        # if model doesn't support nan values, replace with zero to ignore feature
        ###
        if not supports_nan(model):
            X = X.fillna(0)
        
        scores[(tis,ensg)] = pd.Series(model.predict(X), index = X.index)
        missing_count[(tis,ensg)] = {'missing_features':missing_ct,'total_features':len(model.feature_names_in_), 'missing_feature%':missing_ct/len(model.feature_names_in_)}
        
        ###
        # describe feature importances
        # absmax shap value within this cohort
        # save to feature_summary
        ###
        try:
            explainer = shap.TreeExplainer(model, feature_perturbation='tree_path_dependent' )
            shap_values = explainer(X)
        except shap.utils._exceptions.InvalidModelError:
            masker = shap.maskers.Independent(X, max_samples=X.shape[0])
            explainer = shap.Explainer(model, masker)
            shap_values = explainer(X)
        except IndexError:
            assert type(model) is xgb.sklearn.XGBRegressor and model.n_estimators == 0 #null xgb model with zero trees
            shap_values = X.copy()*0
        
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
        # if expression data provided, evaluate prediction accuracy
        ###

        try:
            olap = sorted(set(columns)&set(expr.index))
            ypred = pd.Series(model.predict(X), index = X.index).loc[olap].astype(float)
            y = expr[ensg].loc[olap]
            model_performance[ensg] = {'r2':r2_score(y,ypred), 'rmse':rmse(y,ypred),
                                         'same score frequency':ypred.round(6).value_counts().max()/ypred.shape[0],
                                         'k':X.shape[1], 'n':ypred.shape[0]}
            if ypred.round(6).nunique()>1:
                for corr_method in ['spearman','pearson']:
                    model_performance[ensg][f'{corr_method} corr'] = y.corr(ypred, method = corr_method)
        
            
                ols = sm.OLS(y, sm.add_constant(pd.DataFrame({'EAGLES':ypred}).round(6)))
                ols_results = ols.fit()
                ols_results = ols_results.summary2().tables[0][[2,3]].set_index(2)[3]
            
                model_performance[ensg]['F'] = ols_results['F-statistic:']
                model_performance[ensg]['F_pval'] = ols_results['Prob (F-statistic):']
        except KeyError:
            pass
        
    

    

            
    if args.outdir == "":
        outbase = ""
    else:
        outbase = f'{args.outdir}/'
        
    if len(scores)>0:
        scores = pd.DataFrame(scores)
        scores.to_csv(f'{outbase}scores.tsv', sep = '\t')
        
    if len(missing_count)>0:
        missing_count = pd.DataFrame(missing_count).T
        missing_count.to_csv(f'{outbase}missing_counts.tsv', sep = '\t')
        
    if len(feature_summary)>0:    
        feature_summary = pd.concat(feature_summary).reset_index()
        feature_summary.columns = ['TISSUE','GENE','SNP','SHAP']
        feature_summary.to_csv(f'{outbase}feature_summary.tsv', sep = '\t', index = False)
    
    if len(model_performance)>0:
        model_performance = pd.DataFrame(model_performance).T.astype({'k':int,'n':int})

        for col in ['F','F_pval']:
            model_performance[col] = model_performance[col].replace('nan',np.nan)
        
        r2_adj = (1 - ((1 - model_performance['r2'])*(model_performance['n'] - 1)/(model_performance['n'] - model_performance['k'] - 1).replace(0,np.nan))).sort_values()
        r2_adj_mask = (((model_performance['n'] - model_performance['k'] - 1)<=0)|(r2_adj<0)).replace({True:np.nan, False:1.0})
        model_performance['r2_adj'] = r2_adj*r2_adj_mask
        model_performance.to_csv(f'{outbase}model_performance.tsv', sep = '\t')


if __name__ == "__main__":
    main()
