#!/usr/bin/env python

"""
Date: June 2nd 2021
Author: Alba Sanchis-Juan

This script reformats a BED file to VCF and makes empty genotype 
fields for all the samples provided in the input file
"""

import argparse
import pandas

def main():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('--bed', dest='bed', help='BED input file')
    parser.add_argument('--header', dest='header', help='Header')
    parser.add_argument('--samples', dest='samples', help='Samples')
    parser.add_argument('--vcf', dest='vcf', help='VCF output file')
    args = parser.parse_args()

    bedin = args.bed
    header = args.header
    vcfout = args.vcf
    samples = args.samples

    """
    Get samples
    """
    sampFile=open(samples, 'r')
    sample_ids=sampFile.read().splitlines()
    sample_ids_list = '\t'.join(sample_ids)

    """
    Print header
    """
    out=open(vcfout, 'w')
    with open(header, 'r') as lines:
        for line in lines:
            if line.startswith("##"):
                out.write(line)
            elif line.startswith("#CHROM"):
                ref_line=line.rstrip() + '\t' + sample_ids_list + '\n'
                out.write(ref_line)
    out.close()

    """
    Reformat bed to vcf
    """
    colNames = ['id', 'chrom', 'start', 'end', 'type', 'samples']
    bed = pandas.read_csv(bedin, sep='\t', header=None, names=colNames)

    bed['len'] = bed['end'] - bed['start'] + 1

    bed['ref'] = 'N'
    bed['alt'] = '.'
    bed['qual'] = 'PASS'
    bed['formatVcf'] = 'GT:GQ:RD_CN:RD_GQ:PE_GT:PE_GQ:SR_GT:SR_GQ:EV'

    for samp in sample_ids:
        bed[samp] = './.:.:.:.:.:.:.:.:.'

    bed['typeRef'] = '<' + bed['type'] + '>'

    bed['info'] = 'END=' + bed['end'].astype(str) + \
                ';SVTYPE=' + bed['type'].astype(str) + \
                ';CHR2=' + bed['chrom'].astype(str) + \
                ';SVLEN=' + bed['len'].astype(str) + \
                ';SAMPLES=' + bed['samples']

    colOut = ['chrom', 'start', 'id', 'ref', 'typeRef', 'alt', 'qual', 'info', 'formatVcf'] + sample_ids
    bedOut = bed[colOut]

    bedOut.to_csv(vcfout, mode = 'a', index=False, sep='\t', header=False)

if __name__ == '__main__':
    main()
