#!/usr/bin/env python
import pandas as pd
import joblib
import os
import argparse
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import xgboost as xgb
from pgenlib import PgenReader
import numpy as np


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



def load_covariates(path):
    cov = pd.read_csv(path, sep="\t")
    if "#IID" in cov.columns:
        cov = cov.rename(columns={"#IID": "IID"})
    cov = cov.set_index("IID")
    cov.index = cov.index.astype(str).str.strip()
    return cov

def load_qtl_table(path, gene_id):
    df = pd.read_csv(path, sep="\t")

    df["ENSG_clean"] = df["ENSG"].str.replace(r"\.\d+$", "", regex=True)
    df = df[df["ENSG_clean"] == gene_id]

    df = df[["SNP", "slope"]].set_index("SNP")
    return df

def apply_flipping(X, qtl_df):
    flip_mask = {}

    common_snps = X.columns.intersection(qtl_df.index)
    #Xf = X.copy()

    for snp in common_snps:
        slope = qtl_df.loc[snp, "slope"]
        if slope < 0:
            #Xf[snp] = 2 - Xf[snp]
            flip_mask[snp] = True
        else:
            flip_mask[snp] = False

    #for snp in X.columns:
    #    if snp not in flip_mask:
    #        flip_mask[snp] = False

    #return Xf, flip_mask
    return flip_mask

def clean_gene_id(gene_id):
    gene_base = gene_id.split('.')[0]
    gene_base = gene_base.split('_')[0]
    return gene_base

def load_data(pgen_path, psam_path, pvar_path, expression_path, covar_path, samples = None):
    X = load_pgen_data(pgen_path, psam_path, pvar_path)
    cov = load_covariates(covar_path)
    
    y = pd.read_csv(expression_path, sep="\t", index_col=0)
    if y.shape[1] == 1:
        y = y.iloc[:, 0]
    
    if samples is not None:
        common_samples = list(set(y.index)&set(X.index)&set(cov.index)&samples)
    else:
        common_samples = list(set(y.index)&set(X.index)&set(cov.index))
        
    X = X.loc[common_samples]
    y = y.loc[common_samples]
    cov = cov.loc[common_samples]
    return X, y, cov

def pca_transform(g, thres = 0.999):
    pca = PCA()
    pc_df = pd.DataFrame(pca.fit_transform(g), index = g.index)
    
    pc_comp = pd.DataFrame(pca.components_, columns = g.columns, index = [f'PC{i}' for i in range(1,pca.components_.shape[0]+1)]).T

    pc_df = pc_df.loc[:, range((np.cumsum(pca.explained_variance_ratio_)<thres).sum())]
    pc_df.columns = [f'PC{i}' for i in range(1,pc_df.shape[1]+1)]
    return pc_df, pc_comp.loc[:, pc_df.columns]

def fit_PCR(X_scaled, cov_scaled, y, scaler, thres):
    
    pcs,loadings = pca_transform(pd.DataFrame(X_scaled, index = y.index, columns = scaler.feature_names_in_[:X_scaled.shape[1]]), thres)

    pc_model = LinearRegression()
    pc_model.fit(pcs.join(cov_scaled), y)

    pc_weights = pd.Series(dict(zip(pc_model.feature_names_in_[:pcs.shape[1]], pc_model.coef_[:pcs.shape[1]])))
    temp = pd.Series(dict(zip(cov_scaled.columns, pc_model.coef_[pcs.shape[1]:])))
    temp.index = temp.index.str.replace('cov|','')
    coef_scaled = pd.concat([loadings.dot(pc_weights), temp])

    snp_model = LinearRegression()
    snp_model.coef_ = coef_scaled.values
    snp_model.intercept_ = pc_model.intercept_
    snp_model.feature_names_in_ = scaler.feature_names_in_
    
    return snp_model

def fit_model(X, y, cov, model_type, thres = 1, qtl_df=None, gene_id=None):
    scaler = StandardScaler()
    temp = X.join(cov)
    X_cov_scaled = pd.DataFrame(scaler.fit_transform(temp), index = temp.index, columns = temp.columns)
    feat_list = list(temp.columns)
    del temp

    if model_type == "elasticnet":
        model = ElasticNetCV(l1_ratio=0.5, cv=5, max_iter=10000)
    elif model_type == "ridge":
        model = RidgeCV()
    elif model_type == "rf":
        model = RandomForestRegressor(n_estimators=100)
    elif model_type == "xgb":
        model = xgb.XGBRegressor(objective="reg:squarederror", n_estimators=100)
    elif model_type == "pcr":
        return {"scaler":scaler, "model":fit_PCR(X_cov_scaled.loc[:, X.columns], X_cov_scaled.loc[:, cov.columns].rename({x:'cov|'+x for x in cov.columns},axis = 1), y, scaler, thres), "feature_names":feat_list}
    elif model_type == "flipallele":
        if qtl_df is None:
            raise ValueError("flipallele model requires --qtl input")
        X_snps = X[[c for c in X.columns if c in qtl_df.index]].copy()
        #Xf, flip_mask = apply_flipping(X_snps, qtl_df)
        flip_mask = apply_flipping(X_snps, qtl_df)
        return {
            "feature_names": X_snps.columns.tolist(),
            "flip_mask": flip_mask,
            "gene": gene_id
        }
    else:
        raise ValueError(f"Unsupported model type: {model_type}")

    model.fit(X_cov_scaled, y)
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
    parser.add_argument("--covariates", required=True, help="Path to covariate file")
    parser.add_argument("--model", choices=["elasticnet", "ridge", "rf", "xgb", "pcr", "flipallele"], default="elasticnet", help="Model type")
    parser.add_argument("--output", required=True, help="Output file path to save trained model")
    parser.add_argument("--gene", required=True, help="Gene ID to model")
    parser.add_argument("--samples", required=True, help="training sample list, one per line")
    parser.add_argument("--thres", required=False, help="for pcr regression, limits number of PCs", default = 1)
    parser.add_argument("--qtl", required=False, help="Path to QTL slope file")

    args = parser.parse_args()
    args.gene = clean_gene_id(args.gene)
    
    samples = set(pd.read_csv(args.samples, header = None)[0])
    X, y_all, cov = load_data(args.pgen, args.psam, args.pvar, args.expression, args.covariates, samples)

    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if args.gene not in y_all.columns:
            raise ValueError(f"Gene {args.gene} not found in expression file.")
        y = y_all[args.gene]

    # load QTL table for flipped option
    qtl_df = load_qtl_table(args.qtl, args.gene) if args.model == "flipallele" else None

    model = fit_model(X, y,cov, args.model, args.thres, qtl_df=qtl_df, gene_id=args.gene)
    joblib.dump(model, args.output)
    print(f"Model saved to: {args.output}")

if __name__ == "__main__":
    main()
