#!/usr/bin/env python

import pickle
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pbase", required=True, help="path to base of plink2 files")
    args = parser.parse_args()

    pvar = f'{args.pbase}.pvar'
    snp_index = {}
    
    with open(pvar, 'r') as file:
        h = '##'
        while h.startswith('##'):
            h = next(file)
    
        h = h.split('\t')
        pos_index = h.index('ID')
    
        for i,line in enumerate(file):
            var_id = line.strip().split('\t')[pos_index]
            snp_index[var_id] = i
    
    with open( f'{args.pbase}.pkl', 'wb') as file:
        pickle.dump(snp_index, file)    

    

if __name__ == "__main__":
    main()