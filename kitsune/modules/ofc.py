import argparse as ap
import collections
import sys
from operator import itemgetter

import numpy as np
from tqdm import tqdm

from kitsune.modules import kitsunejf as jf


def read_params(args):
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog="kitsune (ofc)",
        description=(
            "Calculate an observe feature frequency"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--filenames", nargs="+", required=True, help="Genome files in fasta format")
    p.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    p.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    p.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
    p.add_argument("-t", "--thread", type=int, default=1)
    p.add_argument("-o", "--output", type=str, help="Output filename")
    return p.parse_args(args)


def cal_ofc(fsas, kmers, **karg):
    """ Calculate shannon entropy of observed feature occurrences (OFC)
        ofc(l) = -sum(p ln p)
    Args:
        fsas a list of genome files
        kmers a list of kmer to calculate

    Returns: float shannon entropy

    """

    result = dict()
    
    for kmer in kmers:
        keys = list()
        for fsa in fsas:
            keys.extend(list(jf.Kmercount(fsa, kmer, **karg).keys()))
            
        count_feature = list(collections.Counter((collections.Counter(keys).values())).values())
        lnp = np.log2(count_feature) - (kmer * 2)
        result[kmer] = np.sum(np.exp2(lnp) * lnp) * -1

    return result


def run(args):
    """
    Calculate an observe feature frequency
    """

    # Load command line parameters
    args = read_params(args)

    outdata = cal_ofc(args.filenames, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])

    print(outdata, file=open(args.output, "w+") if args.output else None)


if __name__ == "__main__":
    run(sys.argv)
