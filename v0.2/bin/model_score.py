#!/usr/bin/env python

import pandas as pd
import numpy as np
import joblib
import argparse
from fit_model import load_covariates

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load features
    X = pd.read_csv(args.features, sep="\t", index_col=0)

    # Load covariates
    cov = pd.read_csv(args.covariates, sep="\t")
    if "#IID" in cov.columns:
        cov = cov.rename(columns={"#IID": "IID"})
    cov = cov.set_index("IID")
    cov.index = cov.index.astype(str).str.strip()

    # Align features and covariates
    common_samples = X.index.intersection(cov.index)
    X = X.loc[common_samples]
    cov = cov.loc[common_samples]

    # Concatenate covariates
    X = pd.concat([X, cov], axis=1)

    # Load model
    model_dict = joblib.load(args.model)

    scaler = model_dict["scaler"]
    model = model_dict["model"]
    feature_names = model_dict["feature_names"]

    # Ensure all features present
    missing = set(feature_names) - set(X.columns)
    if missing:
        raise ValueError(f"Missing features in input data: {missing}")

    # Order columns
    X_ordered = X[feature_names]

    # Scale and predict
    X_scaled = scaler.transform(X_ordered)
    y_pred = model.predict(X_scaled)

    output_df = pd.DataFrame({
        "sample_id": X_ordered.index,
        "predicted_expression": y_pred
    }).set_index("sample_id")    
    
    output_df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
