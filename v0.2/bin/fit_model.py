import pandas as pd
import joblib
import os
import argparse
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

def load_data(features_path, expression_path, covariates_path):
    X = pd.read_csv(features_path, sep="\t", index_col=0)
    y = pd.read_csv(expression_path, sep="\t", index_col=0)
    c = pd.read_csv(covariates_path, sep="\t", index_col=0)

    
    X = X.loc[y.index]

    if covariates_path:
        covars = pd.read_csv(covariates_path, sep="\t", index_col=0)
        covars = covars.loc[y.index]
        X = pd.concat([X, covars], axis=1)

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
    return (scaler, model)  
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--covariates", default=None)
    parser.add_argument("--model", choices=["elasticnet", "ridge", "rf", "xgb"], default="elasticnet")
    parser.add_argument("--output", required=True)
    parser.add_argument("--gene", required=True)
    
    args = parser.parse_args()
    
    X, y_all = load_data(args.features, args.expression, args.covariates)
    
    if args.gene not in y_all.columns:
        raise ValueError(f"Gene {args.gene} not found in expression file.")
    
    y = y_all[args.gene]
    
    model = fit_model(X, y, args.model)
    joblib.dump(model, args.output)
    print(f"Model saved to: {args.output}")

if __name__ == "__main__":
    main()

#features_path = "genotype_samples_with_expression.csv"
#expression_path = "expression_subjects_matched.csv"
#output_dir = "/cellar/users/nopopko/projects/eagles/model_testing/models"
#model_type = "elasticnet"
