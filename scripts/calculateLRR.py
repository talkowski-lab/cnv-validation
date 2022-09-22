#!/usr/bin/env python

"""
Date: June 2nd 2021
Author: Alba Sanchis-Juan

This script reformats UKBB array data and calculates exp(LRR)
"""

import argparse
import pandas
import numpy as np

def f(x):
    if x == ".":
        return x
    else:
        return np.exp2(float(x))

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--input', dest='input', help='Input file')
    parser.add_argument('--output', dest='output', help='Output file')
    # parser.add_argument('--samples', dest='samples', help='Samples')
    args = parser.parse_args()

    input = args.input
    output = args.output
    # samples = args.samples

    # """
    # Get samples
    # """
    # sampFile=open(samples, 'r')
    # sample_ids=sampFile.read().splitlines()
    # sample_ids_list = '\t'.join(sample_ids)

    """
    Reformat columns
    """
    df = pandas.read_csv(input, sep='\t', compression='gzip')

    df.columns = df.columns.str.replace(pat='\[.*?\]', repl='', n=1).str.replace(pat=':LRR', repl='', n=1).str.replace(pat='# ', repl='', n=1)

    df.insert(2, 'START', df['POS'])
    df.rename(columns={'POS': 'END'}, inplace = True)

    # df = df[df.iloc[:,4] != "."]  # drop missing LRR values

    """
    Apply exp function
    """
    df.iloc[:,4] = df.iloc[:,4].apply(f)

    """
    Write to output
    """
    df.to_csv(output, mode='a', index=False, sep='\t', header=True)

if __name__ == '__main__':
    main()
