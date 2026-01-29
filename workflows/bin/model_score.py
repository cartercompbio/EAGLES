#!/usr/bin/env python

import pandas as pd
import numpy as np
import joblib
import argparse
from fit_model import load_covariates, load_pgen_data, AlleleCount

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--model", required=True)
    parser.add_argument("--covariates", required=False, default = None, help="Path to covariate file")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load features
    X = load_pgen_data(args.pgen, args.psam, args.pvar)

    # Load covariates
    cov = load_covariates(args.covariates)

    # Align features and covariates
    if cov is not None:
        common_samples = X.index.intersection(cov.index)
        X = X.loc[common_samples]
        cov = cov.loc[common_samples]
    
        # Concatenate covariates
        X = pd.concat([X, cov], axis=1)

    # Load model
    model_dict = joblib.load(args.model)

    try:
        scaler = model_dict["scaler"]
    except KeyError:
        scaler = None

    model = model_dict["model"]
    feature_names = model_dict["feature_names"]
    
    missing = set(feature_names) - set(X.columns)
    for col in missing:
        X[col] = np.nan

    # Order columns
    X_ordered = X[feature_names]

    # Scale and predict
    if scaler is not None:
        X_scaled = scaler.transform(X_ordered)
    else:
        X_scaled = X_ordered
    y_pred = model.predict(X_scaled)

    output_df = pd.DataFrame({
        "sample_id": X_ordered.index,
        "predicted_expression": y_pred
    }).set_index("sample_id")    
  
    output_df.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    main()
