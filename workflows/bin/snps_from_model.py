#!/usr/bin/env python

import joblib
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()
    
    model_dict = joblib.load(args.model)
    with open(args.output, 'w') as file:
        file.write('\n'.join(set(model_dict['model'].feature_names_in_) - set(['AGE','ISFEMALE'])))
    