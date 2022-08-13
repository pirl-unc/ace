#!/usr/bin/python3

"""
The purpose of this python3 script is to generate an ELIspot configuration
given three inputs:
1. Number of total peptides.
2. Number of peptides per pool.
3. Total coverage.

Last updated date: August 12, 2022

Author: Jin Seok (Andy) Lee, Dhuvarakesh Karthikeyan
"""


import argparse
from acelib.logger import get_logger
from acelib.solver import *


logger = get_logger(__name__)


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""Generates an ELIspot configuration."""
    )
    arg_parser.add_argument(
        "--num_peptides",
        dest="num_peptides",
        type=int,
        required=True,
        help="Total number of peptides."
    )
    arg_parser.add_argument(
        "--num_peptides_per_pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool. Please make sure this integer is a factor of the total number of peptides."
    )
    arg_parser.add_argument(
        "--num_coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    arg_parser.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    args = arg_parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    df_solutions = solve_assay_configuration(
        n_peptides=args.num_peptides,
        n_peptides_per_pool=args.num_peptides_per_pool,
        n_coverage=args.num_coverage
    )
    df_solutions.to_csv(args.output_tsv_file, sep='\t', index=False)
