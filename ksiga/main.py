#!/usr/bin/env python

""" Main entry.


"""

import argparse
import sys
import os
import gzip
import pathlib

import numpy as np
from sklearn.preprocessing import normalize

from ksiga import logutil
from ksiga import fsig


def openner(filename, **kwargs):
    """Try to return a sensible filehandle

    Args:
        arg1 (TODO): TODO

    Returns: TODO

    """
    if filename.endswith(".gz"):
        return gzip.open(filename, **kwargs)
    else:
        return open(filename, **kwargs)


def main():
    commands = {"index": index,
                "relent": relative_entropy,
                "cre_kmer": cumulative_relative_entropy,
                "acf": average_common_feature,
                "acf_kmer": acf_kmer,
                "ofc": observe_feature_frequency,
                "dmatrix": generate_distance_matrix,
                "uniq": generate_uniq_mer}

    parser = argparse.ArgumentParser(description="Signature for virus",
                                     usage="""ksiga <command> [<args>]

Commands can be:
index <filenames>                     Compute k-mer.
relent <filename.sig>                 Compute relative entropy.
cre <filename.sig>                    Compute cumulative relative entropy.
acf <filenames.sig>                   Compute average number of common feature between signatures.
ofc <filenames.sig>                   Compute observed feature frequencies.
dmatrix <filenames.sig>               Compute distance matrix.
""")
    parser.add_argument('command')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands:
        parser.print_help()
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
    parser.add_argument("-o", "--wd", default=os.getcwd())
    args = parser.parse_args(args)

    filenames = args.filenames
    ksize = args.ksize
    wd = args.wd

    for filename in args.filenames:
        if not os.path.exists(filename):
            # TODO: Warn or exit here.
            pass

    if not os.path.isdir(wd):  #  Check for output folder
        try:
            os.makedirs(wd)
        except FileExistsError as err:
            logutil.notify("Filename is already exists.")
            exit(1)

    for filename in filenames:
        # Clean folder name from file
        basename = pathlib.Path(filename).name
        outputName = "{wd}/{fn}.ksig".format(wd=wd, fn=basename)
        fInputH = openner(filename, mode="rt")
        fsig.build_signature(fInputH, ksize, outputName)


def relative_entropy(args):
    """ Calculate relative entropy of genome.

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


def cumulative_relative_entropy(args):
    """ Calculate optimal k-mer through CRE value.

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate cumulative relative entropy of all genomes"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-ks", "--kfrom", required=True, type=int, help="Calculate from k-mer")
    parser.add_argument("-ke", "--kend", required=True, type=int, help="last k-mer")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(args)

    filenames = args.filenames
    kmerStart = args.kfrom
    kmerEnd = args.kend

    cre, kmer = fsig.calculate_cre(filenames[0], kmerStart, kmerEnd)
    print(kmer)


def average_common_feature(args):
    """ Calculate an average number of common feature pairwise
        between one genome against others

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate average number of common feature"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-o", "--output")
    args = parser.parse_args(args)

    outF = args.output
    if outF is None:
        outHandle = sys.stdout.buffer
    else:
        outHandle = open(outF, "wb")  # wb for numpy write

    acf = fsig.calculate_average_common_feature(args.filenames, args.ksize)

    np.savetxt(outHandle, acf)


def acf_kmer(args):
    """ Calculate an average number of common feature pairwise
        between one genome against others

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate average number of common feature"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-ks", "--kfrom", required=True, type=int, help="Calculate from k-mer")
    parser.add_argument("-ke", "--kend", required=True, type=int, help="last k-mer")
    parser.add_argument("-o", "--output")
    args = parser.parse_args(args)

    filenames = args.filenames
    outF = args.output
    kmerStart = args.kfrom
    kmerEnd = args.kend

    if outF is None:
        outHandle = sys.stdout.buffer
    else:
        outHandle = open(outF, "wb")  # wb for numpy write

    acf = fsig.calculate_acf(filenames, kmerStart, kmerEnd)

    # np.savetxt(outHandle, acf)


def observe_feature_frequency(args):
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-w", "--wd", default=os.getcwd())
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args(args)

    ksize = args.ksize
    output = args.output

    prob, occ, kmerStr = fsig.calculate_obsff(args.filenames, ksize)
    df = pd.DataFrame({"kmer": kmerStr, "prob": prob, "occ": occ})
    df.set_index("kmer", inplace=True)
    df.to_csv(output, sep="\t")

def generate_uniq_mer(args):
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-w", "--wd", default=os.getcwd())
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args(args)

    ksize = args.ksize
    output = args.output

    allPossible, uniq = fsig.calculate_uniq_mer(args.filenames, ksize)

    outF = args.output
    if outF is None:
        outHandle = sys.stdout.buffer
    else:
        outHandle = open(outF, "wb")  # wb for numpy write

    np.savetxt(outHandle, [allPossible, uniq], fmt="%i")


def generate_distance_matrix(args):
    """Generate distance matrix base on k-mer

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    import ksiga.fsig as fsig
    from ksiga import distance

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-o", "--output")
    parser.add_argument("-t", "--n_thread", type=int, default=1)
    parser.add_argument("-d", "--distance", default="euclid")
    args = parser.parse_args(args)

    fn = distance.DISTANCE_FUNCTION.get(args.distance, None)

    if fn is None:
        allowDistance = list(distance.DISTANCE_FUNCTION.keys())
        print(USAGE)
        exit(1)

    filenames = args.filenames
    ksize = args.ksize
    outF = args.output
    n_thread = args.n_thread

    if outF is None:
        outHandle = sys.stdout.buffer
    else:
        outHandle = open(outF, "wb")  # wb for numpy write

    # Check for existence of file.
    for filename in args.filenames:
        if not os.path.exists(filename):
            # TODO: Do something about this
            pass

    csr_matrix = fsig.rebuild_sparse_matrix(filenames, ksize)
    rowNum = csr_matrix.shape[0]

    csr_matrix_norm = normalize(csr_matrix, norm='l1', axis=1)

    result = fn(csr_matrix_norm)
    np.savetxt(outHandle, result)
    logutil.notify("Result is written to {}".format(outF))
    sys.exit(0)
