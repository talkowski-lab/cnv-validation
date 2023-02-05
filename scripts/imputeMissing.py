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

    df = pd.read_csv(input, sep='\t', na_values=['.'], index_col='ID', compression='gzip')

    sample_cols = df.columns[3:]

    null_samples = df.iloc[:, 3:].isna()

    val_order = df.iloc[:, 3:].values.argsort()

    probes_with_missing_vals = null_samples.any(axis=1)

    with open(missing, 'w') as missing_report:
        for probe in probes_with_missing_vals.index[probes_with_missing_vals]:
            probe_nas = null_samples.loc[probe]
            null_samples_for_row = sample_cols[probe_nas]
            missing_report.write("\t".join(
                map(str,
                    df.loc[probe, ['CHROM', 'START', 'END']].tolist() +
                    [probe] +
                    [",".join([sample for sample in null_samples_for_row]) + "\n"])))
            notnull = len(sample_cols) - null_samples.loc[probe].sum()
            for null_sample in null_samples_for_row:
                x = rng.integers(0, notnull)
                f = sample_cols[val_order[df.index.get_loc(probe)] == x]
                newval = df.loc[probe, f].item()
                df.loc[probe, null_sample] = newval

    df.to_csv(output, mode='w', index=True, sep='\t', header=True, compression='gzip')


if __name__ == '__main__':
    main()
