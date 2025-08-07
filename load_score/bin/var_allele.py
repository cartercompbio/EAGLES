#!/usr/bin/env python

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="specify the variant allele for plink counting")
    parser.add_argument("--fpath", required=True, help="Path to ldpruned snps")

    args = parser.parse_args()

    with open(args.fpath, 'r') as infile:
        snps = infile.read().strip().split('\n')
        print('\n'.join(map(lambda x: x+'\t'+x[x.rfind(':')+1:], snps)))
