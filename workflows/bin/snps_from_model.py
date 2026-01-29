#!/usr/bin/env python

import joblib
import argparse
from fit_model import AlleleCount
from qtl_filter import variant_ids_from_pvar

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--covariates", required=False, default=None)
    parser.add_argument("--pvar", required=True)
    
    args = parser.parse_args()
    covariate_names = set()
    if args.covariates is not None:
        with open(args.covariates, 'r') as file:
            l = next(file).strip()
            covariate_names = set(l.split('\t')) | set(l.split(','))
        
    model_dict = joblib.load(args.model)
    model_snps = set(model_dict['model'].feature_names_in_) - covariate_names
     
    if len(model_snps)>0:
        with open(args.output, 'w') as file:
            file.write('\n'.join(model_snps))