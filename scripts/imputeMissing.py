#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--input', dest='input', help='Input tsv')
    parser.add_argument('--output', dest='output', help='Output tsv')
    parser.add_argument('--missing-report', dest='missing', help='BED file containing missing probe data locations')
    parser.add_argument('--seed', dest='seed', help='Random Seed', default=42, type=int)
    args = parser.parse_args()

    input = args.input
    output = args.output
    missing = args.missing
    seed = args.seed

    rng = np.random.default_rng(seed=seed)

    df = pd.read_csv(input, sep='\t', compression='gzip', na_values=['.'], index_col='ID')

    sample_cols = df.columns[3:]

    with open(missing, 'w') as missing_report:
        for probe in df.index:
            null_samples = df.columns[df.loc[probe].isnull()]
            if len(null_samples) > 0:
                missing_report.write("\t".join(
                    map(str,
                        df.loc[probe, ['CHROM', 'START', 'END']].tolist() +
                        [probe] +
                        [",".join([sample for sample in df.columns[df.loc[probe].isnull()]]) + "\n"])))
                notnull = df.shape[1] - 3 - len(null_samples)
                val_order = df.loc[probe][3:].argsort()
                for null_sample in null_samples:
                    x = rng.integers(0, notnull)
                    f = sample_cols[val_order == x]
                    newval = df.loc[probe, f].item()
                    df.loc[probe, null_sample] = newval

    df.to_csv(output, mode='w', index=False, sep='\t', header=True)


if __name__ == '__main__':
    main()
