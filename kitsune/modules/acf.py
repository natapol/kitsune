import argparse as ap
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
        prog="kitsune (acf)",
        description=(
            "Calculate an average number of common feature pairwise between one genome against others"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--filenames", nargs="+", type=str, required=True, help="Genome files in fasta format"),
    p.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    p.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    p.add_argument("-k", "--kmers", nargs="+", required=True, type=int, help="Have to state before")
    p.add_argument("-t", "--thread", type=int, default=1)
    p.add_argument("-o", "--output", type=str, help="Output filename")
    return p.parse_args(args)


def cal_acf(fsas, kmers, **karg):
    """Calculate Average number of common features (ACF)

    Args:
        fsas (str):  genome file name(s).
        kmers (): a list of kmer to calculate.

    Kwargs:
        state (bool): Current state to be in.
        thread (int): Number of thread to calculate default 1
        lower (int): default 1
        bchashsize (str): hashsize for jellyfish bc step default '1G'
        hashsize (str): hashsize for jellyfish count step default '100M'
        canonical (bool): set canonical calculation

    Returns:
        dict(kmer: acf)

    Raises:
        AttributeError, KeyError

    A really great idea.  A way you might use me is

    >>> print public_fn_with_googley_docstring(name='foo', state=None)
    0

    BTW, this always returns 0.  **NEVER** use with :class:`MyPublicClass`.

    """
    
    n = len(fsas)
    
    if n >= 2:
        result = dict()
        
        for kmer in tqdm(kmers):
            keys_array = list()
            
            for fsa in fsas:
                keys_array.append(set(jf.Kmercount(fsa, kmer, **karg).keys()))
            
            ccf = 0
            
            for pri_idx, key in enumerate(keys_array):
                for sec_idx in range(pri_idx + 1, n):
                    ccf += len(keys_array[pri_idx] & keys_array[sec_idx])
            
            result[kmer] = ccf/(n-1)
        
        return result
    
    else:
        raise


def run(args):
    """
    Calculate an average number of common feature pairwise
    between one genome against others
    """

    # Load command line parameters
    args = read_params(args)

    outdata = cal_acf(args.filenames, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])

    print(outdata, file=open(args.output, "w+") if args.output else None)


if __name__ == "__main__":
    run(sys.argv)
