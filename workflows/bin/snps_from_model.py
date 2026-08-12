#!/usr/bin/env python

from fit_model import AlleleCount
import joblib
import argparse

def get_feats(model_path):
    model_dict = joblib.load(model_path)
    return list(model_dict['model'].feature_names_in_)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, help="path to .pkl model file")
    parser.add_argument("--outfile", required=True, help="path to write snps")

    args = parser.parse_args()

    feats = get_feats(args.model)
    with open(args.outfile, 'w') as file:
        file.write('\n'.join(feats))