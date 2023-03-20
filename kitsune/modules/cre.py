import argparse as ap
import math
import sys
from operator import itemgetter

from tqdm import tqdm

from kitsune.modules import kitsunejf as jf


def read_params(args):
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog="kitsune (cre)",
        description=(
            "Calculate k-mer from cumulative relative entropy of all genomes"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--filename", type=str, required=True, help="A genome file in fasta format")
    p.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    p.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    p.add_argument("-ke", "--kend", required=True, type=int, help="Last k-mer")
    p.add_argument("-kf", "--kfrom", default=4, type=int, help="Calculate from k-mer")
    p.add_argument("-t", "--thread", type=int, default=1)
    p.add_argument("-o", "--output", type=str, help="Output filename")
    return p.parse_args(args)


def cal_cre(fsa, kend, kfrom, **karg):
    """ Calculate Cumulative Relative Entropy (CRE)
        CRE = sum(RE from kmer to infinite)
    Args:
        fsa genome file in fasta format
        kend the infinite kmer
        kfrom calculate to (defualt=4)

    Returns: dict(kmer: CRE)

    """
    
    a0 = None
    a1 = None
    a2 = None
    
    result = dict()
    
    for kmer in tqdm(range(kend, kfrom-1, -1)):
        if a0 is None:
            a0 = jf.Kmercount(fsa, kmer, **karg)
            a1 = jf.Kmercount(fsa, kmer - 1, **karg)
            a2 = jf.Kmercount(fsa, kmer - 2, **karg)
        
        else:
            a0 = a1
            a1 = a2
            a2 = jf.Kmercount(fsa, kmer - 2, **karg)
        
        if kmer + 1 in result:
            result[kmer] = cal_re(a0, a1, a2) + result[kmer + 1]
        
        else:
            result[kmer] = cal_re(a0, a1, a2)
    
    return result


def cal_re(a0, a1, a2):
    """ Calculate Relative Entropy (RE)
        f' = fl * fr / fb
        example: f' of mer 'ATTTGCA' 
                f  is frequency of ATTTGCA from a0
                fl is frequency of ATTTGC- from a1
                fr is frequency of -TTTGCA from a1
                fb is frequency of -TTTGC- from a2
        RE = f * log2(f/f')
    Args:
        a0 Kmercount kmer
        a1 Kmercount kmer - 1
        a2 Kmercount kmer - 2

    Returns: float

    """

    result = 0
    rfactor = a0.sum
    lfactor = math.log((a1.sum ** 2) / (a2.sum * rfactor), 2)
    
    for key in a0.keys():
        realf = a0[key]
        left = a1[key[:-1]] if key[:-1] in a1 else 0
        right = a1[key[1:]] if key[1:] in a1 else 0
        below = a2[key[1:-1]] if key[1:-1] in a2 else 0
        
        if 0 not in (left, right, below):
            expectf = left * right / below
            result += max(0, realf / rfactor * (math.log(realf / expectf, 2) + lfactor))
    
    return result


def run(args):
    """
    Calculate k-mer from cumulative relative entropy of all genomes
    """

    # Load command line parameters
    args = read_params(args)

    outdata = cal_cre(args.filename, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])

    print(outdata, file=open(args.output, "w+") if args.output else None)


if __name__ == "__main__":
    run(sys.argv)
