#!/usr/bin/env python


import argparse
import sys
import os
import pathlib

from ksiga import fsig

DEFAULT_K = 13

def main():
    commands = {"index": index, "relent": relative_entropy, "acf":average_common_feature, "off": off}
    parser = argparse.ArgumentParser(description="Signature for virus",
                                     usage="""ksiga <command> [<args>]

Commands can be:
index <filenames>                      Compute k-mer.
relent <filenames.sig>                 Compute relative entropy.
acf <filename1.sig> <filename2.sig>    Compute average number of common feature between signature1 and signature2.
off <filenames.sig>                    Compute observed feature frequencies.
""")
    parser.add_argument('command')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands:
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])

def index(args):
    """ Create index for input sequences

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="+", help="file(s) of sequences")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-w", "--wd", default=os.getcwd())
    args = parser.parse_args(args)

    for filename in args.filenames:
        # Clean folder name from file
        basename = pathlib.Path(filename).name
        os.path.exists(filename)
        outputName = "{wd}/{fn}.hdf".format(wd=args.wd, fn=basename)
        fsig.build_signature(filename, args.ksize, outputName)

def relative_entropy(args):
    """ Calculate relative entropy of

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--file", required=True, help="")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-w", "--wd", default=os.getcwd())
    args = parser.parse_args(args)

    relEntropy = fsig.calculate_relative_entropy(args.file, args.ksize)
    print(relEntropy)


def average_common_feature(args):
    """ Calculate an average number of common feature.
        aka, k-mer exist in both sample.

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in", required=True, help="")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    exit(1)

def off(args):
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    pass

