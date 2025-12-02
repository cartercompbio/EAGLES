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
    Xf = X.copy()

    for snp in common_snps:
        slope = qtl_df.loc[snp, "slope"]
        if slope < 0:
            Xf[snp] = 2 - Xf[snp]
            flip_mask[snp] = True
        else:
            flip_mask[snp] = False

    for snp in X.columns:
        if snp not in flip_mask:
            flip_mask[snp] = False

    return Xf, flip_mask

def clean_gene_id(gene_id):
    gene_base = gene_id.split('.')[0]
    gene_base = gene_base.split('_')[0]
    return gene_base

def load_data(features_path, expression_path, samples = None):
    X = pd.read_csv(features_path, sep="\t", index_col=0)
    
    y = pd.read_csv(expression_path, sep="\t", index_col=0)
    if y.shape[1] == 1:
        y = y.iloc[:, 0]
        
    common_samples = y.index.intersection(X.index)
    if samples is not None:
        common_samples = list(set(common_samples)&samples)
        
    X = X.loc[common_samples]
    y = y.loc[common_samples]

    return X, y

def pca_transform(g, thres = 0.999):
    pca = PCA()
    pc_df = pd.DataFrame(pca.fit_transform(g), index = g.index)
    
    pc_comp = pd.DataFrame(pca.components_, columns = g.columns, index = [f'PC{i}' for i in range(1,pca.components_.shape[0]+1)]).T

    pc_df = pc_df.loc[:, range((np.cumsum(pca.explained_variance_ratio_)<thres).sum())]
    pc_df.columns = [f'PC{i}' for i in range(1,pc_df.shape[1]+1)]
    return pc_df, pc_comp.loc[:, pc_df.columns]

def fit_PCR(X_scaled, y, scaler, thres):
    
    pcs,loadings = pca_transform(pd.DataFrame(X_scaled, index = y.index, columns = scaler.feature_names_in_), thres)

    pc_model = LinearRegression()
    pc_model.fit(pcs, y)

    pc_weights = pd.Series(dict(zip(pc_model.feature_names_in_, pc_model.coef_)))
    coef_scaled = loadings.dot(pc_weights)

    coef_orig = coef_scaled / scaler.scale_
    intercept_orig = (pc_model.intercept_ - 
                         np.sum(coef_scaled * scaler.mean_ / scaler.scale_))

    snp_model = LinearRegression()
    snp_model.coef_ = coef_orig
    snp_model.intercept_ = intercept_orig
    snp_model.feature_names_in_ = scaler.feature_names_in_
    
    return snp_model

def fit_model(X, y, model_type, thres = 1, qtl_df=None, gene_id=None):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    if model_type == "elasticnet":
        model = ElasticNetCV(l1_ratio=0.5, cv=5, max_iter=10000)
    elif model_type == "ridge":
        model = RidgeCV()
    elif model_type == "rf":
        model = RandomForestRegressor(n_estimators=100)
    elif model_type == "xgb":
        model = xgb.XGBRegressor(objective="reg:squarederror", n_estimators=100)
    elif model_type == "pcr":
        return {"scaler":scaler, "model":fit_PCR(X_scaled, y, scaler, thres), "feature_names":X.columns.tolist()}
    elif model_type == "flipallele":
        if qtl_df is None:
            raise ValueError("flipallele model requires --qtl input")
        X_snps = X[[c for c in X.columns if c in qtl_df.index]].copy()
        Xf, flip_mask = apply_flipping(X_snps, qtl_df)
        return {
            "feature_names": X_snps.columns.tolist(),
            "flip_mask": flip_mask,
            "gene": gene_id
        }
    else:
        raise ValueError(f"Unsupported model type: {model_type}")

    model.fit(X_scaled, y)
    return {
        "scaler": scaler,
        "model": model,
        "feature_names": X.columns.tolist()
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True, help="Path to features file (e.g. genotype data)")
    parser.add_argument("--expression", required=True, help="Path to expression data file")
    parser.add_argument("--model", choices=["elasticnet", "ridge", "rf", "xgb", "pcr", "flipallele"], default="elasticnet", help="Model type")
    parser.add_argument("--output", required=True, help="Output file path to save trained model")
    parser.add_argument("--gene", required=True, help="Gene ID to model")
    parser.add_argument("--samples", required=True, help="training sample list, one per line")
    parser.add_argument("--thres", required=False, help="for pcr regression, limits number of PCs", default = 1)
    parser.add_argument("--qtl", required=False, help="Path to QTL slope file")

    args = parser.parse_args()
    args.gene = clean_gene_id(args.gene)
    
    samples = set(pd.read_csv(args.samples, header = None)[0])
    X, y_all = load_data(args.features, args.expression, samples)

    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if args.gene not in y_all.columns:
            raise ValueError(f"Gene {args.gene} not found in expression file.")
        y = y_all[args.gene]

    # load QTL table for flipped option
    qtl_df = load_qtl_table(args.qtl, args.gene) if args.model == "flipallele" else None

    model = fit_model(X, y, args.model, args.thres, qtl_df=qtl_df, gene_id=args.gene)
    joblib.dump(model, args.output)
    print(f"Model saved to: {args.output}")

if __name__ == "__main__":
    main()
