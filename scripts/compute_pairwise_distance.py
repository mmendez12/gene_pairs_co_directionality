#!/usr/bin/env python3

import argparse
import pandas as pd

def compute_distances(df):
    res = []
    for chrom, data in df.groupby('#contig'):
        data = data.reset_index(drop=True)
        for i, row1 in data.iterrows():
            for j in range(i + 1, len(data)):
                row2 = data.iloc[j]
                if row1['stop'] < row2['start']:
                    # row2 is downstream of row1
                    distance = row2['start'] - row1['stop']
                elif row1['start'] > row2['stop']:
                    # row2 is upstream of row1
                    distance = row1['start'] - row2['stop']
                else:
                    # intervals are overlapping
                    continue

                res.append({
                    'chrom': chrom,
                    'gene1': row1['name'],
                    'gene2': row2['name'],
                    'distance': distance
                })
    return pd.DataFrame(res)


def main(input_file, output):
    df = pd.read_table(input_file, usecols=['#contig', 'name', 'start', 'stop'])
    result_df = compute_distances(df)
    result_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute pairwise gene distances for non-overlapping genes.")
    parser.add_argument("input", help="Input TSV file (gzipped ok)")
    parser.add_argument("output", help="Output TSV file")

    args = parser.parse_args()
    main(args.input, args.output)
