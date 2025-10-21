import pandas as pd
import joblib
import os
import argparse
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

def clean_gene_id(gene_id):
    gene_base = gene_id.split('.')[0]
    gene_base = gene_base.split('_')[0]
    return gene_base

def load_data(features_path, expression_path):
    X = pd.read_csv(features_path, sep="\t", index_col=0)
    
    y = pd.read_csv(expression_path, sep="\t", index_col=0)
    if y.shape[1] == 1:
        y = y.iloc[:, 0]

    common_samples = y.index.intersection(X.index)
    X = X.loc[common_samples]
    y = y.loc[common_samples]

    return X, y

def fit_model(X, y, model_type="elasticnet"):
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
    parser.add_argument("--model", choices=["elasticnet", "ridge", "rf", "xgb"], default="elasticnet", help="Model type")
    parser.add_argument("--output", required=True, help="Output file path to save trained model")
    parser.add_argument("--gene", required=True, help="Gene ID to model")

    args = parser.parse_args()
    args.gene = clean_gene_id(args.gene)

    X, y_all = load_data(args.features, args.expression)

    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if args.gene not in y_all.columns:
            raise ValueError(f"Gene {args.gene} not found in expression file.")
        y = y_all[args.gene]

    model = fit_model(X, y, args.model)
    joblib.dump(model, args.output)
    print(f"Model saved to: {args.output}")

if __name__ == "__main__":
    main()
