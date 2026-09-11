#!/usr/bin/env python

import pandas as pd
import pickle
import os
import argparse
from functools import reduce

def load_by_id(filename, index_file, gene, col_name = 'SNP'):
    if filename.endswith('.csv'):
        sep = ','
    elif filename.endswith('.tsv'):
        sep = '\t'
    else:
        raise ValueError(f'filename must be .csv or .tsv, not {filename}')
        
    with open(index_file, 'rb') as idx:
        id_ranges = pickle.load(idx)
        
    if gene not in id_ranges:
        return []
        
    start_byte,end_byte = id_ranges[gene]

    with open(filename, 'rb') as infile:
        l = next(infile).decode().strip().split(sep)
        row_pos = l.index(col_name)

        infile.seek(start_byte)

        values = []
        while infile.tell() < end_byte:
            line = infile.readline().decode().strip()
            if line:
                values.append(line.split(sep)[row_pos])

        return values
    
def variant_ids_from_pvar(pvar):
    with open(pvar, 'r') as file:
        lines = [l.strip().replace("#","").split('\t') for l in file if not l.startswith("##")]
        if len(lines)==0:
            return set()
        id_pos = lines[0].index('ID')
        pvar_snps = set([l[id_pos] for l in lines[1:]])
        return pvar_snps
    
def load_by_tis(gene, outfile, **kwargs):
    try:
        qtl_file = kwargs['qtl_file']
    except:
        raise ValueError('qtl_file kwarg missing')
        
    try:
        index_file = kwargs['index_file']
    except:
        raise ValueError('index_folder kwarg missing')
        
    try:
        col_name = kwargs['col_name']
    except:
        col_name = 'SNP'
        
    try:
        pvar_path = kwargs['pvar']
        assert pvar_path.endswith('.pvar')
    except:
        pvar_path = None

    snp_set = set(load_by_id(qtl_file, index_file, gene))
    
    if pvar_path is not None:
        pvar_snps = set(variant_ids_from_pvar(pvar_path))
        merged_set = snp_set&pvar_snps
    
    
    if len(merged_set)>0:
        with open(outfile, 'w') as outfile:
            outfile.write('\n'.join(merged_set))        
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get gene-associated eQTLs")
    parser.add_argument("--qtl-file", required=True, help="Path to qtl table")
    parser.add_argument("--index-file", required=True, help="Path to qtl table index file")
    parser.add_argument("--gene", required=True, help="gene identifier for desired eQTLs")
    parser.add_argument("--tis", required=True, help="tissue")
    parser.add_argument("--output", required=True, help="Path to output eqtl list")
    parser.add_argument("--pvar", default=None, help="Path to .pvar file to restrict returned eQTLs")

    args = parser.parse_args()

    # Load data
    load_by_tis(args.gene, args.output, qtl_file = args.qtl_file, index_file = args.index_file, tis = args.tis, pvar=args.pvar)