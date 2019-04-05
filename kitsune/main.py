#!/usr/bin/env python

""" Main entry of the script.
"""

import argparse
import sys
import os
import gzip
import datetime
import numpy as np
from operator import itemgetter, attrgetter

# import internal library
from . import kitsunejf as jf

USAGE = """
"""

def main():
    commands = {
                "cre": cumulative_relative_entropy,
                "acf": average_common_feature,
                "ofc": observe_feature_occurrence,
                "dmatrix": generate_distance_matrix
               }

    parser = argparse.ArgumentParser(description="Signature for virus",
                                     usage="""ksiga <command> [<args>]

Commands can be:
cre <filename>                    Compute cumulative relative entropy.
acf <filenames>                   Compute average number of common feature between signatures.
ofc <filenames>                   Compute observed feature frequencies.
dmatrix <filenames>               Compute distance matrix.
""")
    parser.add_argument('command')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands:
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])

def cumulative_relative_entropy(args):
    from . import cre
    """ Calculate optimal k-mer through CRE value.

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate k-mer from cumulative relative entropy of all genomes"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", type=str, help="a genome file in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-ke", "--kend", required=True, type=int, help="last k-mer")
    parser.add_argument("-kf", "--kfrom", default=4, type=int, help="Calculate from k-mer")
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)
    
    outdata = cre.cal_cre(args.filename, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)


def average_common_feature(args):
    from . import acf
    """ Calculate an average number of common feature pairwise
        between one genome against others

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate average number of common feature"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    outdata = cal_acf(args.filenames, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

def observe_feature_occurrence(args):
    from . import ofc
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """

    desc = "Calculate observe feature occurrence"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    outdata = cal_ofc(args.filenames, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

def generate_distance_matrix(args):
    from . import matrix
    """Generate distance matrix base on k-mer

    The output will 

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate a distance matrix"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmer", required=True, type=int)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("-d", "--distance", default="cosine")
    parser.add_argument("-f", "--format", default="phylip")
    args = parser.parse_args(args)
    args.dist_func = jf.DISTANCE_FUNCTION[args.distance]
    outdata = matrix.cal_matrix(args.filenames, **vars(args))
    if args.format == "phylip":
        outdata = str(len(outdata)) + '\n' + '\n'.join([k + '\t' + '\t'.join([str(x) for x in v]) for k, v in outdata.items()])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

if __name__ == "__main__":
    main()