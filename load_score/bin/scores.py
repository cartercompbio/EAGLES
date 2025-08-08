#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
from scipy import sparse
from scipy.sparse import coo_matrix
from sklearn.decomposition import PCA

def load_data(G_path, E_path):
    """
    Read in genotype matrix (G) and effect size matrix (E)
    """
    G = pd.read_csv(G_path, sep='\t', index_col=0)
    E = pd.read_csv(E_path, sep='\t', index_col=0)
    return G, E

def matrix_mult_scoring(G, E):
    """
    Calculate gene scores using matrix multiplication.
    
    Inputs:
    G: Genotype matrix (samples x SNPs)
    E: Effect size matrix (SNPs x genes)
    
    Outputs:
    S: Gene scores matrix (samples x genes)
    """
    
    S = G.values @ E.values

    sample_names = G.index
    gene_names = E.columns
    
    S_df = pd.DataFrame(S, columns=gene_names, index=sample_names)
    return S_df

def directional_scoring(G, E):
    """
    Calculate gene scores accounting for directional effects.
    
    Inputs:
    G: Genotype matrix (samples x SNPs)
    E: Effect size matrix (SNPs x genes)
    
    Outputs:
    S: Gene scores matrix (samples x genes)
    """

    sample_names = G.index
    snp_names = G.columns
    gene_names = E.columns

    E_coo = coo_matrix(E.values)

    scores = np.zeros((len(sample_names), len(gene_names)))

    for snp_idx, gene_idx, effect in zip(E_coo.row, E_coo.col, E_coo.data):
        snp_name = E.index[snp_idx]
        if snp_name not in G.columns:
            continue

        alt_counts = G[snp_name].values
        ref_counts = 2 - alt_counts

        if effect > 0:
            contribution = effect * alt_counts
        else:
            contribution = abs(effect) * ref_counts
            # contribution = (1 / abs(effect)) * ref_counts is another possibility

        scores[:, gene_idx] += contribution

    S_df = pd.DataFrame(scores, index=sample_names, columns=gene_names)
    return S_df

def direction_allele_count(G, E):
    """
    Normalized count of allele associated with higher expression
    
    Inputs:
    G: Genotype matrix (samples x SNPs)
    E: Effect size matrix (SNPs x genes)
        
    Outputs:
    S: Gene scores matrix (samples x genes)
        values are between [0,1]
    """
    
    res = {}

    for gene in E:

        pos_snps = E[E[gene]>0].index
        neg_snps = E[E[gene]<0].index

        score = (G[pos_snps].join(2 - G[neg_snps]).sum(axis = 1)/(2*(len(pos_snps) + len(neg_snps))))
        res[gene] = score
    res = pd.DataFrame(res).loc[G.index, E.columns].fillna(0)
    return res

def pc1_score(G, E):
    """
    PC1 computed over all eQTLs associated with a gene
        (count of allele associated with higher expression)
    
    Inputs:
    G: Genotype matrix (samples x SNPs)
    E: Effect size matrix (SNPs x genes)
        
    Outputs:
    S: Gene scores matrix (samples x genes)
        values are between [0,1]
    """
    res = {}

    for gene in (E!=0).any().replace(False, np.nan).dropna().index:

        pos_snps = E[E[gene]>0].index
        neg_snps = E[E[gene]<0].index
        cur = G[pos_snps].join(2 - G[neg_snps])
        try:
            assert cur.shape[1]>1
            pca = PCA(n_components = 1)
            pca.fit(cur.T)
            score = pd.Series(pca.components_[0], index = cur.index)
        except AssertionError:
            score = cur.iloc[:, 0]
        score = (score - score.min())/(score.max() - score.min())

        res[gene] = score
    res = pd.DataFrame(res)
    return res

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Gene scoring")
    parser.add_argument("--genotypes", required=True, help="Path to genotype matrix")
    parser.add_argument("--effects", required=True, help="Path to SNP-gene effect matrix")
    parser.add_argument("--method", choices=["matrixmult", "directional"], default="directional", help="Scoring method: 'matrixmult', 'directional', TBA...")
    parser.add_argument("--output", required=True, help="Path to output gene scores file")

    args = parser.parse_args()

    # Load data
    G_df, E_df = load_data(args.genotypes, args.effects)

    print(f"Genotype matrix shape: {G_df.shape} (samples x SNPs)")
    print(f"Effect matrix shape:   {E_df.shape} (SNPs x genes)")

    # Ensure snps are aligned in same order
    common_snps = G_df.columns.intersection(E_df.index)
    G_df = G_df[common_snps]
    E_df = E_df.loc[common_snps]
    
    print(f"Using {len(common_snps)} SNPs for scoring.")

    # Then proceed with scoring method selected
    if args.method == "matrixmult":
        result = matrix_mult_scoring(G_df, E_df)
    else:
        result = directional_scoring(G_df, E_df)

    result.to_csv(args.output, sep="\t")
    print(f"Gene scores written to: {args.output}")
