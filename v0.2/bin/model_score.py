#!/usr/bin/env python

import pandas as pd
import numpy as np
import joblib
import argparse
from fit_model import load_covariates, load_pgen_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True, help="path to .pgen file")
    parser.add_argument("--psam", required=True, help="path to .psam file")
    parser.add_argument("--pvar", required=True, help="path to .pvar file")
    parser.add_argument("--model", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load features
    X = load_pgen_data(args.pgen, args.psam, args.pvar)

    # Load covariates
    cov = load_covariates(args.covariates)

    # Align features and covariates
    common_samples = X.index.intersection(cov.index)
    X = X.loc[common_samples]
    cov = cov.loc[common_samples]

    # Concatenate covariates
    X = pd.concat([X, cov], axis=1)

    # Load model
    model_dict = joblib.load(args.model)

    # Handle flipped allele model
    if "flip_mask" in model_dict:
        # Only keep SNPs (ignore covariates)
        snp_cols = model_dict["feature_names"]
        X_snps = X[snp_cols].copy() 

        # Apply flips according to flip_mask
        for snp, flip in model_dict["flip_mask"].items():
            if flip:
                X_snps[snp] = 2 - X_snps[snp]

        # Sum alleles per sample
        y_pred = X_snps.sum(axis=1)

        output_df = pd.DataFrame({
            "sample_id": X_snps.index,
            "flipped_allele_sum": y_pred
        }).set_index("sample_id")

    # ML models
    else:
        scaler = model_dict["scaler"]
        model = model_dict["model"]
        feature_names = model_dict["feature_names"]

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
