#!/usr/bin/env python
import pandas as pd
import numpy as np
import joblib
import argparse

# import functions from fit_model.py and model_score.py
from fit_model import (
    clean_gene_id,
    load_data,
    load_qtl_table,
    fit_model,
)

from model_score import (
    load_covariates,
    load_pgen_data,
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgen", required=True)
    parser.add_argument("--psam", required=True)
    parser.add_argument("--pvar", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--model", choices=["elasticnet", "ridge", "rf", "xgb", "pcr", "flipallele"], required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--thres", default=1, help="PCR variance threshold")
    parser.add_argument("--qtl", required=False)
    parser.add_argument("--model_output", required=True)
    parser.add_argument("--score_output", required=True)

    args = parser.parse_args()
    gene_id = clean_gene_id(args.gene)

    samples = set(pd.read_csv(args.samples, header=None)[0])

    X, y_all, cov = load_data(
        args.pgen,
        args.psam,
        args.pvar,
        args.expression,
        args.covariates,
        samples
    )
    if isinstance(y_all, pd.Series):
        y = y_all
    else:
        if gene_id not in y_all.columns:
            raise ValueError(f"{gene_id} not found in expression file.")
        y = y_all[gene_id]

    # Load qtl table for flipped option
    qtl_df = load_qtl_table(args.qtl, gene_id) if args.model == "flipallele" else None

    # Fit model process
    model_dict = fit_model(
        X,
        y,
        cov,
        args.model,
        float(args.thres),
        qtl_df=qtl_df,
        gene_id=gene_id
    )
    
    # Save model 
    model_path = f"{args.model_output}.pkl"
    joblib.dump(model_dict, model_path)
    
    # Scoring process
    X_pred = load_pgen_data(args.pgen, args.psam, args.pvar)
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
        }).set_index("sample_id")

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
        }).set_index("sample_id")

    score_path = f"{args.score_output}_scores.tsv"
    output.to_csv(score_path, sep="\t")    
    print("Done: trained and scored", gene_id)

if __name__ == "__main__":
    main()
