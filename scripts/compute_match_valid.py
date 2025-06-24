#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd


def match(v1, v2):
    return (v1 == v2).sum()


def valid(v1, v2):
    return (~np.isnan(v1) & ~np.isnan(v2)).sum()


def main(input_file, match_out, valid_out):
    df = pd.read_table(input_file, index_col=0, low_memory=False)
    match_df = df.T.corr(method=match)
    valid_df = df.T.corr(method=valid)

    match_df.to_csv(match_out, sep='\t')
    valid_df.to_csv(valid_out, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute match and valid correlation matrices")
    parser.add_argument("input")
    parser.add_argument("match")
    parser.add_argument("valid")
    args = parser.parse_args()

    main(args.input, args.match, args.valid)
