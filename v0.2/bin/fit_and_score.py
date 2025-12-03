import pandas as pd
import numpy as np
import joblib
import argparse
from sklearn.linear_model import ElasticNetCV, RidgeCV, LinearRegression
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

def clean_gene_id(gene_id):
    gene_base = gene_id.split('.')[0]
    gene_base = gene_base.split('_')[0]
    return gene_base

def load_data(features_path, expression_path, samples=None):
    X = pd.read_csv(features_path, sep="\t", index_col=0)
    y = pd.read_csv(expression_path, sep="\t", index_col=0)

    if y.shape[1] == 1:
        y = y.iloc[:, 0]

    common_samples = X.index.intersection(y.index)

    if samples is not None:
        common_samples = list(set(common_samples) & samples)

    X = X.loc[common_samples]
    y = y.loc[common_samples]
    return X, y

def load_qtl_table(path, gene_id):
    if path is None:
        return None
    df = pd.read_csv(path, sep="\t")
    df["ENSG_clean"] = df["ENSG"].str.replace(r"\.\d+$", "", regex=True)
    df = df[df["ENSG_clean"] == gene_id]
    return df[["SNP", "slope"]].set_index("SNP")

def apply_flipping(X, qtl_df):
    flip_mask = {}
    common_snps = X.columns.intersection(qtl_df.index)
    for snp in common_snps:
        slope = qtl_df.loc[snp, "slope"]
        flip_mask[snp] = (slope < 0)
    return flip_mask

def pca_transform(X, thres=0.999):
    pca = PCA()
    pcs = pd.DataFrame(pca.fit_transform(X), index=X.index)

    loadings = pd.DataFrame(
        pca.components_,
        columns=X.columns,
        index=[f"PC{i}" for i in range(1, pca.components_.shape[0] + 1)]
    ).T

    keep = (np.cumsum(pca.explained_variance_ratio_) < thres).sum()
    pcs = pcs.iloc[:, :keep]
    pcs.columns = [f"PC{i}" for i in range(1, pcs.shape[1] + 1)]
    return pcs, loadings.loc[:, pcs.columns]

# Model fitting
def fit_model(X, y, model_type, thres=1, qtl_df=None, gene_id=None):
    # Flipped allele scenario
    if model_type == "flipallele":
        X_snps = X[[c for c in X.columns if c in qtl_df.index]].copy()
        flip_mask = apply_flipping(X_snps, qtl_df)
        return {
            "feature_names": X_snps.columns.tolist(),
            "flip_mask": flip_mask,
            "gene": gene_id
        }

    # Standard ML model scenario
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
        pcs, loadings = pca_transform(pd.DataFrame(X_scaled, index=y.index, columns=X.columns), thres)
        pc_model = LinearRegression().fit(pcs, y)

        pc_weights = pd.Series(pc_model.coef_, index=pcs.columns)
        coef_scaled = loadings.dot(pc_weights)
        coef_orig = coef_scaled / scaler.scale_
        intercept_orig = pc_model.intercept_ - np.sum(coef_scaled * scaler.mean_ / scaler.scale_)

        snp_model = LinearRegression()
        snp_model.coef_ = coef_orig
        snp_model.intercept_ = intercept_orig
        snp_model.feature_names_in_ = X.columns

        return {
            "scaler": scaler,
            "model": snp_model,
            "feature_names": X.columns.tolist()
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
    parser.add_argument("--features", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--qtl", required=False)
    parser.add_argument("--thres", default=1)
    parser.add_argument("--out-prefix", required=True)

    args = parser.parse_args()
    gene_id = clean_gene_id(args.gene)


    train_ids = set(pd.read_csv(args.samples, header=None)[0])

    X, y_all = load_data(args.features, args.expression, train_ids)
    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if gene_id not in y_all.columns:
            raise ValueError(f"{gene_id} not found in expression!")
        y = y_all[gene_id]

    qtl_df = load_qtl_table(args.qtl, gene_id) if args.model == "flipallele" else None

    # Fit model
    model_dict = fit_model(X, y, args.model, args.thres, qtl_df=qtl_df, gene_id=gene_id)

    # Save model 
    joblib.dump(model_dict, f"{args.out_prefix}.pkl")

    # Scoring with same data
    X_pred = pd.read_csv(args.features, sep="\t", index_col=0)
    cov = load_covariates(args.covariates)
    
    common = X_pred.index.intersection(cov.index)
    X_pred = X_pred.loc[common]
    cov = cov.loc[common]
    X_pred = pd.concat([X_pred, cov], axis=1)

    # Flipped allele scenario
    if "flip_mask" in model_dict:
        snps = model_dict["feature_names"]
        X_snps = X_pred[snps].copy()

        for snp, flip in model_dict["flip_mask"].items():
            if flip:
                X_snps[snp] = 2 - X_snps[snp]

        y_pred = X_snps.sum(axis=1)

        output = pd.DataFrame({
            "sample_id": X_snps.index,
            "flipped_allele_sum": y_pred
        })

    # Normal ML scenaio
    else:
        scaler = model_dict["scaler"]
        model = model_dict["model"]
        feats = model_dict["feature_names"]

        X_ordered = X_pred[feats]
        X_scaled = scaler.transform(X_ordered)
        y_pred = model.predict(X_scaled)

        output = pd.DataFrame({
            "sample_id": X_ordered.index,
            "predicted_expression": y_pred
        })

    output.set_index("sample_id").to_csv(f"{args.out_prefix}_scores.tsv", sep="\t")
    print("Done: trained and scored", gene_id)

if __name__ == "__main__":
    main()
