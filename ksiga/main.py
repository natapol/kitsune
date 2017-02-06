#!/usr/bin/env python


import argparse
import sys
import os
import gzip
import pathlib

from ksiga import fsig

DEFAULT_K = 13

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
                "acf":average_common_feature,
                "off": observe_feature_frequency,
                "dmatrix": generate_distance_matrix}

    parser = argparse.ArgumentParser(description="Signature for virus",
                                     usage="""ksiga <command> [<args>]

Commands can be:
index <filenames>                     Compute k-mer.
relent <filename.sig>                 Compute relative entropy.
acf <filenames.sig>                   Compute average number of common feature between signatures.
off <filenames.sig>                   Compute observed feature frequencies.
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
            # TODO: Warn?
            pass

    for filename in filenames:
        # Clean folder name from file
        basename = pathlib.Path(filename).name
        outputName = "{wd}/{fn}.ksig".format(wd=wd, fn=basename)
        fInputH = openner(filename, mode="rt")
        fsig.build_signature(fInputH, ksize, outputName)


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
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    args = parser.parse_args(args)
    acf = fsig.calculate_average_common_feature(args.filenames, args.ksize)
    print(acf)


def observe_feature_frequency(args):
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--file", required=True, help="")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-w", "--wd", default=os.getcwd())
    args = parser.parse_args(args)
    occF = fsig.calculate_obsff(args.file, args.ksize)
    print(occF)


def generate_distance_matrix(args):
    """Generate distance matrix base on k-mer

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    import ksiga.fsig as fsig
    import ksiga.mmath as mmath
    import itertools
    from ksiga.ksignature import KmerSignature

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="+", help="file(s) of signature")
    parser.add_argument("-k", "--ksize", required=True, type=int)
    parser.add_argument("-o", "--output")
    args = parser.parse_args(args)

    filesname = args.filenames
    ksize = args.ksize
    outF = args.output

    if outF is None:
        outHandle = sys.stdout
    else:
        outHandle = open(outF, "w")

    # Check for existence of file.
    for filename in args.filenames:
        if not os.path.exists(filename):
            # TODO: Do something about this
            pass

    csr_matrix = fsig.rebuild_sparse_matrix(args.filenames, args.ksize)
    rowNum = csr_matrix.shape[0]

    # TODO: Maybe use pairwise_distances from scikit-learn?
    for i, j in itertools.combinations(range(rowNum), r=2):
        iRow = csr_matrix[i]
        jRow = csr_matrix[j]
        distance = mmath.sparse_js_distance(iRow, jRow)
        outHandle.write(str(distance))
        outHandle.write(os.linesep)
    print("Finish")
