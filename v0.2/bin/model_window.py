#!/usr/bin/env python

import pandas as pd
import argparse

def get_gene_info(info_path, upstream, downstream, **kwargs):
    '''
    upstream/downstream: parameters defining regulatory window around genes
    
    info_path is to a csv/tsv containing information about gene locations
    must contain at least 4 columns: 
        1. gene_name (default "ENSG")
        2. chrom (default "CHROM")
        3. strand (default "STRAND")
        4. gene_start (default "TSS")
        
    if something other than default column name is used, must be specified as kwarg (e.g. gene_name = "HGNC")
    
    
    returns new dataframe including above 4 columns plus:
        1. "start"
        2. "stop"
        
    where chrom:start-stop is a region including the gene and its regulatory region
    '''
    try:
        gene_name = kwargs['gene_name']
    except KeyError:
        gene_name = 'ENSG'
        
    try:
        chrom = kwargs['chrom']
    except:
        chrom = 'CHROM'
        
    try:
        strand = kwargs['strand']
    except:
        strand = 'STRAND'
    
    try:
        gene_start = kwargs['gene_start']
    except KeyError:
        gene_start = 'TSS'
    
    if info_path.endswith(','):
        sep = ','
    else:
        sep = '\t'
    
    df = pd.read_csv(info_path, sep = sep)[[gene_name, chrom, strand, gene_start]]

    df_plus = df[df[strand]=='+'].copy()
    df_neg = df[df[strand]=='-'].copy()
    del df


    for df_sign,up,down in zip([df_plus,df_neg],[upstream, downstream],[downstream,upstream]):
        temp = (df_sign[gene_start] - up).copy()
        temp[temp<0] = 0
        df_sign['start'] = temp
        df_sign['end'] = df_sign[gene_start] + down
        
    return pd.concat([df_plus,df_neg]).sort_index()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True)
    parser.add_argument("--upstream", type=int, required=True)
    parser.add_argument("--downstream", type=int, required=True)
    parser.add_argument("--output", required=True)


    args = parser.parse_args()
    gene_info = get_gene_info(args.path, args.upstream, args.downstream)
    gene_info.to_csv(f'{args.output}.tsv', sep = '\t', index = False)
