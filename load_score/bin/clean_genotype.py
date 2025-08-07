#!/usr/bin/env python

import pandas as pd
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="clean plink .raw file into a more simple .tsv format")
    parser.add_argument("--inpath", required=True, help="Path to plink .raw genotype file")
    parser.add_argument("--outpath", required=True, help="Path to .tsv genotype output file")
    
    args = parser.parse_args()
    
    df = pd.read_csv(args.inpath, sep = '\t', index_col = 'IID').drop(['FID','PAT','MAT','SEX','PHENOTYPE'],axis = 1)
    df.columns = df.columns.str.split('_',expand = True).get_level_values(0)

    df.to_csv(args.outpath, sep = '\t')
