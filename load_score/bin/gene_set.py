import pandas as pd
import argparse

def load_long_effect(E_path):
    """
    Load long-form SNP-gene-effect table with columns: Gene, Variant, slope (could allow for flexibility later on)
    """
    df = pd.read_csv(E_path, sep='\t')
    if not {'Gene', 'Variant', 'slope'}.issubset(df.columns):
        raise ValueError("Input SNP x Gene file must contain 'Gene', 'Variant', and 'slope' columns.")
    return df

def load_gene_set(gene_set_path):
    """
    Load gene set from file (one gene per line)
    """
    with open(gene_set_path) as f:
        return [line.strip() for line in f if line.strip()]

def to_dense_matrix(df_long, gene_set):
    """
    Convert long-form SNP-gene-effect table to dense matrix, keeping only genes in the gene set.
    """
    df_filtered = df_long[df_long["Gene"].isin(gene_set)]
    df_dense = df_filtered.pivot(index="Variant", columns="Gene", values="slope")
    return df_dense.fillna(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create dense SNP x gene table only including genes from input gene set")
    parser.add_argument("--effects", required=True, help="Path to SNP-gene-effect table")
    parser.add_argument("--gene-set", required=True, help="Path to gene set file")
    parser.add_argument("--output", required=True, help="Path to output")

    args = parser.parse_args()

    df_long = load_long_effect(args.effects)
    gene_set = load_gene_set(args.gene_set)

    df_dense = to_dense_matrix(df_long, gene_set)
    df_dense.to_csv(args.output, sep='\t')
    print(f"Filtered SNP × Gene matrix saved to: {args.output}")
