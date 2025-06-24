#!/usr/bin/env python3

import pandas as pd
import argparse

def main(chrom, distances_file, match_file, valid_file, output_file):
    # Load distance matrix
    distance = pd.read_table(distances_file)
    distance = distance[distance.chrom == chrom].copy()

    for tag, path in (("match", match_file), ("valid", valid_file)):
        df = pd.read_table(path, index_col=0)
        flat_df = df.stack(future_stack=True).reset_index()
        flat_df = flat_df.rename(columns={
            'name': 'gene1',
            'level_1': 'gene2',
            0: tag
        }).dropna()

        distance = distance.merge(flat_df, on=['gene1', 'gene2'], how='left')
        distance.dropna(inplace=True)

    distance.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge distances with match/valid data.")
    parser.add_argument('--chrom', required=True, help="Chromosome name (e.g. chr1)")
    parser.add_argument('--distances', required=True, help="Path to pairwise distances file")
    parser.add_argument('--match', required=True, help="Path to match score file")
    parser.add_argument('--valid', required=True, help="Path to valid score file")
    parser.add_argument('--output', required=True, help="Output path for merged TSV")

    args = parser.parse_args()
    main(args.chrom, args.distances, args.match, args.valid, args.output)