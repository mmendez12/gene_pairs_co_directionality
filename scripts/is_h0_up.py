#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.stats import binomtest
import argparse

def is_h0_up(element):
    if pd.isna(element):
        return element
    h0, h1 = map(int, element.split('|'))
    total = h0 + h1
    if total <= 8:
        return np.nan
    elif binomtest(h0, total).pvalue >= 0.05:
        return np.nan
    return h0 > h1

def main(input_file, output_file):
    df = pd.read_table(input_file, index_col='name')
    df = df[[col for col in df.columns if col.startswith('GTEX')]]
    df = df.replace('0|0', np.nan)
    result = df.map(is_h0_up)
    result.to_csv(output_file, sep='\t', )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine if haplotype 0 is upregulated")
    parser.add_argument("input", help="Input phased genotype file (gzipped)")
    parser.add_argument("output", help="Output CSV file")
    args = parser.parse_args()
    main(args.input, args.output)