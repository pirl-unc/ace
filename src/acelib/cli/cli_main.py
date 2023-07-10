# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement the primary ACE command.
"""


import argparse
import acelib
from .cli_generate import *
from .cli_deconvolve import *
from .cli_verify import *


logger = get_logger(__name__)


def init_arg_parser():
    """
    Initializes the input argument parser.

    Returns
    -------
    An instance of argparse.ArgumentParser
    An instance of argparse.ArgumentParser subparsers
    """
    arg_parser = argparse.ArgumentParser(
        description="ACE Configurator for ELISpot."
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(acelib.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='ACE sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_ace_generate_arg_parser(sub_parsers=sub_parsers)      # generate
    sub_parsers = add_ace_deconvolve_arg_parser(sub_parsers=sub_parsers)    # deconvolve
    sub_parsers = add_ace_verify_arg_parser(sub_parsers=sub_parsers)        # verify
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'generate':
        run_ace_generate_from_parsed_args(args=args)
    elif args.which == 'deconvolve':
        run_ace_deconvolve_from_parsed_args(args=args)
    elif args.which == 'verify':
        run_ace_verify_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)

