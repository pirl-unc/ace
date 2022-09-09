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
from acelib.solver import *
from acelib.visualization import *


logger = get_logger(__name__)


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""Generates an ELIspot configuration."""
    )
    arg_parser._action_groups.pop()
    arg_parser_required = arg_parser.add_argument_group('required arguments')
    arg_parser_required.add_argument(
        "--num_peptides",
        dest="num_peptides",
        type=int,
        required=True,
        help="Total number of peptides."
    )
    arg_parser_required.add_argument(
        "--num_peptides_per_pool",
        dest="num_peptides_per_pool",
        type=int,
        required=True,
        help="Number of peptides per pool. Please make sure this integer is a factor of the total number of peptides."
    )
    arg_parser_required.add_argument(
        "--num_coverage",
        dest="num_coverage",
        type=int,
        required=True,
        help="Total coverage (i.e. number of peptide replicates)."
    )
    arg_parser_required.add_argument(
        "--num_threads",
        dest="num_threads",
        type=int,
        required=True,
        default=2,
        help="Number of threads to parallelize the computation (default: 2). Recommended: 8."
    )
    arg_parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    arg_parser_optional = arg_parser.add_argument_group('optional arguments')
    arg_parser_optional.add_argument(
        "--output_pdf_file",
        dest="output_pdf_file",
        type=str,
        required=False,
        help="Output PDF file."
    )
    args = arg_parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    # Step 1. Generate assay configuration
    df_configuration = generate_assay_configuration(
        n_peptides=args.num_peptides,
        n_peptides_per_pool=args.num_peptides_per_pool,
        n_coverage=args.num_coverage,
        num_threads=args.num_threads
    )

    # Step 2. Save to TSV file
    df_configuration.to_csv(args.output_tsv_file, sep='\t', index=False)

    # Step 3. Save to PDF file
    if args.output_pdf_file is not None:
        plot_configuration_table(df_configuration=df_configuration,
                                 save_figure=True,
                                 output_pdf_file=args.output_pdf_file)
