#!/usr/bin/env python

import argparse
import joblib
import pandas as pd
import numpy as np
from pgenlib import PgenReader
import pickle

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression, ElasticNetCV
from sklearn.decomposition import PCA

import optuna
import xgboost as xgb



def load_slopes(filename, index_file, gene, **kwargs):
    '''
    kwargs:
        var_col: string, which column from filename table to use as variant index
        val_col: string, which column from filename table to use as eqtl slope values
        
    returns pandas.Series
    
    '''
    
    try:
        var_col = kwargs['var_col']
    except KeyError:
        var_col = 'SNP'
        
    try:
        val_col = kwargs['val_col']
    except:
        val_col = 'slope'

    if filename.endswith('.csv'):
        sep = ','
    elif filename.endswith('.tsv'):
        sep = '\t'
    else:
        raise ValueError(f'filename must be .csv or .tsv, not {filename}')

    with open(index_file, 'rb') as idx:
        id_ranges = pickle.load(idx)

    if gene not in id_ranges:
        return pd.Series()

    start_byte,end_byte = id_ranges[gene]

    with open(filename, 'rb') as infile:
        l = next(infile).decode().strip().split(sep)
        infile.seek(start_byte)

        values = []
        while infile.tell() < end_byte:
            line = infile.readline().decode().strip()
            if line:
                values.append(line.split(sep))
                
        return pd.DataFrame(values, columns = l).set_index(var_col)[val_col]


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
        
    genotype_matrix = pd.DataFrame(genotype_matrix, index = variants, columns = samples).T
        
    return genotype_matrix

def clean_gene_id(gene_id):
    gene_base = gene_id.split('.')[0]
    gene_base = gene_base.split('_')[0]
    return gene_base

def load_data(pgen_path, psam_path, pvar_path, expression_path, samples = None, **kwargs):
    try:
        gene=kwargs['gene']
    except:
        raise ValueError('missing kwarg "gene"')
    
    X = load_pgen_data(pgen_path, psam_path, pvar_path)

    with open(expression_path, 'r') as file:
        gene_col_index = next(file).strip().split('\t').index(gene)
    
    y = pd.read_csv(expression_path, sep="\t", index_col=0, usecols=[0,gene_col_index])
    if y.shape[1] == 1:
        y = y.iloc[:, 0]
    
    common_samples = list(set(y.index)&set(X.index)&set(samples))
        
    y = y.loc[common_samples]
    X = X.loc[y.index]

    return X, y

def objective(trial, X, y, model_type, penalty=0.1, n_splits=5, early_stopping_rounds=50):
    if model_type != 'xgb':
        raise ValueError('valid options are "xgb"')

    params = {
        'objective': 'reg:squarederror',
        'eval_metric': 'rmse',
        'booster': 'gbtree',
        'tree_method': 'hist',  # required for grow_policy='lossguide'

        'lambda': trial.suggest_float('lambda', 0.1, 1000, log=True),
        'alpha': trial.suggest_float('alpha', 0.1, 1000, log=True),

        'max_depth': trial.suggest_int('max_depth', 2, 10),
        'eta': trial.suggest_float('eta', 0.01, 0.5, log=True),
        'gamma': trial.suggest_float('gamma', 1e-2, 1000, log=True),
        'grow_policy': trial.suggest_categorical('grow_policy', ['depthwise', 'lossguide']),

        'subsample': trial.suggest_float('subsample', 0.2, 0.8),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 0.8),
        'colsample_bylevel': trial.suggest_float('colsample_bylevel', 0.5, 0.8),

        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'n_estimators': 2000, 
        'random_state': 100,
    }

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=100)
    train_scores = []
    val_scores = []
    best_iters = []

    X_arr = np.asarray(X)
    y_arr = np.asarray(y)

    for train_idx, val_idx in kf.split(X_arr):
        X_tr, X_val = X_arr[train_idx], X_arr[val_idx]
        y_tr, y_val = y_arr[train_idx], y_arr[val_idx]

        model = xgb.XGBRegressor(
            **params,
            early_stopping_rounds=early_stopping_rounds,
        )
        model.fit(
            X_tr, y_tr,
            eval_set=[(X_val, y_val)],
            verbose=False,
        )

        best_iters.append(model.best_iteration)

        train_pred = model.predict(X_tr)
        val_pred = model.predict(X_val)

        train_scores.append(np.sqrt(np.mean((y_tr - train_pred) ** 2)))
        val_scores.append(np.sqrt(np.mean((y_val - val_pred) ** 2)))

    train_score = np.mean(train_scores)
    val_score = np.mean(val_scores)
    overfit_gap = max(val_score - train_score, 0)

    # stash the actual tree counts so you can inspect/reuse them later
    trial.set_user_attr('best_iterations', best_iters)
    trial.set_user_attr('mean_best_iteration', float(np.mean(best_iters)))

    return val_score + penalty * overfit_gap


def tune_model(X, y, model_type):
    study = optuna.create_study(
        direction='minimize',
        sampler=optuna.samplers.TPESampler(seed=100),
    )
    study.optimize(lambda trial: objective(trial, X, y, model_type),
                   n_trials=200, show_progress_bar=False)

    best_trial = study.best_trial
    best_params = dict(best_trial.params)
    best_params['n_estimators'] = round(best_trial.user_attrs['mean_best_iteration'])

    return best_params
    
def fit_flipallele(X, eqtl):
    model = LinearRegression()
    qtl_ser = eqtl[eqtl.index.isin(X.columns)].astype(float)
    
    model.intercept_ = int(2*(qtl_ser<0).sum())
    model.coef_ = np.array(2*(qtl_ser>0) - 1).astype(int)
    model.feature_names_in_ = np.array(qtl_ser.index)

    return model

def pca_transform(X_scaled, thres = 0.999, **kwargs):
    trial_number = 1
    while trial_number < 5:
        try:
            pca = PCA()
            pc_df = pd.DataFrame(pca.fit_transform(X_scaled), index = X_scaled.index)
            pc_comp = pd.DataFrame(pca.components_, columns = X_scaled.columns, index = [f'PC{i}' for i in range(1,pca.components_.shape[0]+1)]).T
            
            n_comps = min([1 + list((1-np.cumsum(pca.explained_variance_ratio_)).round(6) > (1-thres)).index(False),
                           pc_df.shape[0]])
            
            pc_df = pc_df[range(n_comps)]
            pc_df.columns = [f'PC{i}' for i in range(1,pc_df.shape[1]+1)]
            pc_comp = pc_comp.loc[:, pc_df.columns]
            return pc_df, pc_comp

        except:
            trial_number += 1

            

def fit_PCR(X_scaled, y, scaler, thres):
    
    pcs,loadings = pca_transform(X_scaled, thres)
    
    # Check if PCA returned any components
    if pcs.shape[1] == 0:
        return None

    pc_model = LinearRegression()
    pc_model.fit(pcs, y)
    pc_weights = pd.Series(dict(zip(pc_model.feature_names_in_[:pcs.shape[1]], pc_model.coef_[:pcs.shape[1]])))
    coef_scaled = loadings.dot(pc_weights)

    snp_model = LinearRegression()
    snp_model.coef_ = coef_scaled.values
    snp_model.intercept_ = pc_model.intercept_
    snp_model.feature_names_in_ = scaler.feature_names_in_
    
    return snp_model   

def fit_model(X, y, model_type, thres = 1, qtl_ser=None, gene_id=None):
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), index = X.index, columns = X.columns)
    feat_list = list(X.columns)
    
    if model_type == "flipallele":
        model = fit_flipallele(X, qtl_ser)
        return {"model":model, "feature_names":feat_list}
    elif X_scaled.shape[1]==1:
        model = LinearRegression()
    elif model_type == "elasticnet":
        model = ElasticNetCV(l1_ratio=0.5, cv=5, max_iter=10000)
    elif model_type == "xgb":
        best_params = tune_model(X_scaled, y, model_type)
        model = xgb.XGBRegressor(**best_params, random_state = 100)
    elif model_type == "pcr":
        model = fit_PCR(X_scaled, y, scaler, thres)
            
        if model is not None:
            return {"scaler":scaler, "model":model, "feature_names":feat_list}
        else:
            return None

    else:
        raise ValueError(f"Unsupported model type: {model_type}")

    model.fit(X_scaled, y)
    return {
        "scaler": scaler,
        "model": model,
        "feature_names": feat_list
    }



def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument("--features", required=True, help="Path to features file (e.g. genotype data)")
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--expression", required=True, help="Path to expression data file")
    parser.add_argument("--model", choices=["elasticnet", "xgb", "pcr", "flipallele"], default="elasticnet", help="Model type")
    parser.add_argument("--output", required=True, help="Output file path to save trained model")
    parser.add_argument("--gene", required=True, help="Gene ID to model")
    parser.add_argument("--samples", required=True, help="training sample list, one per line")
    parser.add_argument("--thres", required=False, help="for pcr regression, limits number of PCs", default = 1.0, type=float)
    parser.add_argument("--qtl", required=False, help="Path to QTL slope file")
    parser.add_argument("--qtl-index", required=False, help="Path to QTL slope index file")

    args = parser.parse_args()
    args.gene = clean_gene_id(args.gene)
    
    samples = set(pd.read_csv(args.samples, header = None)[0])
    X, y_all = load_data(args.pgen, args.psam, args.pvar, args.expression, samples, gene=args.gene)

    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if args.gene not in y_all.columns:
            raise ValueError(f"Gene {args.gene} not found in expression file.")
        y = y_all[args.gene]

    # load QTL table for flipped option
    qtl_ser = load_slopes(args.qtl, args.qtl_index, args.gene).astype(float) if args.model == "flipallele" else None
    if qtl_ser is not None:
        qtl_ser = qtl_ser[qtl_ser.index.isin(X)]

    model = fit_model(X, y, args.model, args.thres, qtl_ser=qtl_ser, gene_id=args.gene)
    if model is not None:
        joblib.dump(model, args.output)
        print(f"Model saved to: {args.output}")

if __name__ == "__main__":
    main()
