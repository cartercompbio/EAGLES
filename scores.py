import numpy as np
import pandas as pd
import argparse
from scipy import sparse
from scipy.sparse import coo_matrix

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Gene scoring")
    parser.add_argument("--genotypes", required=True, help="Path to genotype matrix")
    parser.add_argument("--effects", required=True, help="Path to SNP-gene effect matrix")
    parser.add_argument("--method", choices=["matrixmult", "directional"], default="directional", help="Scoring method: 'matrixmult', 'directional', TBA...")
    parser.add_argument("--output", required=True, help="Path to output gene scores file")

    args = parser.parse_args()

    # Load data
    G_df, E_df = load_data(args.genotypes, args.effects)

    if args.method == "matrixmult":
        result = matrix_mult_scoring(G_df, E_df)
    else:
        result = directional_scoring(G_df, E_df)

    result.to_csv(args.output, sep="\t")
    print(f"Gene scores written to: {args.output}")
