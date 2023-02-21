import argparse as ap
import collections
import math
import multiprocessing as mp
import os
import time

from collections import Counter
from itertools import combinations, chain, islice
from operator import itemgetter, attrgetter

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
            "Find the recommended optimal kmer using acf,cre and ofc for a given set of genomes. "
            "Example: kitsune kopt genomeList.txt -kl 15 --canonical --fast -t 4 -o out.txt"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--filenames",
        type=valid_file,
        required=True,
        help="A file that list the path to all genomes (fasta format) with extension as (.txt, .csv, .tab) or no extension"
    )
    p.add_argument("--fast", action='store_true', help="Jellyfish one-pass calculation (faster)")
    p.add_argument("--canonical", action='store_true', help="Jellyfish count only canonical mer")
    p.add_argument(
        "-kl",
        "--klarge",
        type=int,
        required=True,
        help="Largest k-mer length to consider, note: the smallest kmer length is 4"
    )
    p.add_argument("-o", "--output", type=str, help="Output filename")
    p.add_argument("--closely_related", action='store_true', help="For closely related set of genomes, use this option" )
    p.add_argument(
        "-x",
        "--cre_cutoff",
        type=number,
        default=0.1,
        help="Cutoff to use in selecting kmers whose cre's are <= (cutoff * max(cre)), Default = 10 percent, ie x=0.1"
    )
    p.add_argument(
        "-y",
        "--acf_cutoff",
        type=number,
        default=0.1,
        help="Cutoff to use in selecting kmers whose acf's are >= (cutoff * max(acf)), Default = 10 percent, ie y=0.1"
    )
    p.add_argument("-t", "--thread", type=int, default=1, help="Number of threads")
    p.add_argument("-n", "--nproc", type=int, default=1, help="Number of processes")
    return p.parse_args(args)


def valid_file(filepath):
    _, ext = os.path.splitext(filepath)
    
    if ext.lower() not in ('.csv', '.tab', '.txt', ''):
        raise argparse.ArgumentTypeError('File must have a csv, txt, tab extension or no extension')
    
    return filepath


def number(value):
    try:
        value = float(value)

    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (value,))

    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (value,))

    return value


def cal_cre(fsa, kl, **karg):
    """ Calculate Cumulative Relative Entropy (CRE)
        CRE = sum(RE from kmer to infinite)
    Args:
        fsa genome file in fasta format
        kend the infinite kmer
        kfrom calculate to (defualt=4)

    Returns: dict(kmer: CRE)

    """
    ksmall = 4
    klarge = kl
    a0 = None
    a1 = None
    a2 = None
    result = {}
    for kmer in tqdm(range(klarge,ksmall-1,-1)):
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
        left = a1[key[0:-1]]
        right = a1[key[1:]]
        below = a2[key[1:-1]]
        if 0 not in (left, right, below):
            expectf = left * right / below
            result += max(0, realf / rfactor * (math.log(realf / expectf, 2) + lfactor))
    return result


def cal_acf(fsas, klarge, **karg):
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
    ksmall = 4
    n = len(fsas)
    if n >= 2:
        result = {}
        for kmer in tqdm(range(ksmall, klarge+1)):
            if kmer % 2 == 1:
                keys_array = []
                for fsa in fsas:
                    keys_array.append(set(jf.Kmercount(fsa, kmer, **karg).keys()))
                ccf = 0
                for pri_idx, key in enumerate(keys_array):
                    for sec_idx in range(pri_idx + 1, n):
                        ccf += len(keys_array[pri_idx] & keys_array[sec_idx])
                result[kmer] = ccf/(n-1)
            else:
                continue
        return result
    else:
        raise


def cal_ofc(fsas, kmer, **karg):
    """ Calculate shannon entropy of observed feature occurrences (OFC)
        ofc(l) = -sum(p ln p)
    Args:
        fsa a genome file
        k_mers a list of kmer to calculate

    Returns: float shannon entropy

    """
    result = {}

    keys = []
    if kmer % 2 == 1:
        for fsa in fsas:
            keys.extend(list(jf.Kmercount(fsa, kmer, **karg).keys()))

        count_feature = list(collections.Counter((collections.Counter(keys).values())).values())
        lnp = np.log2(count_feature) - (2 * kmer)
        result[kmer] = np.sum(np.exp2(lnp) * lnp) * -1
    return result


def filter(n, m): 
    c = Counter(n)
    cr = set(n)
    result = []
    
    for i in cr:
        if c[i] >= int(m * 0.5):
            result.append(i)
    
    return result


def kopt_cre(l, kl, genomeList, args):
    genome = genomeList[l]
    outdata = cal_cre(genome, kl, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    return outdata


def kopt_acf(p,args): 
    outdata_acf = cal_acf(p, **vars(args))
    outdata_acf = sorted(outdata_acf.items(), key=itemgetter(0))
    return outdata_acf


def kopt_ofc(genomeList,kmer, args):
    outdata_ofc = cal_ofc(genomeList,kmer, **vars(args))
    outdata_ofc = sorted(outdata_ofc.items(), key=itemgetter(0))
    return outdata_ofc


def run(args):
    """
    Find the recommended optimal kmer using acf,cre and ofc for a given set of genomes
    """
    
    # Load command line parameters
    args = read_params(args)

    genomeList = open(args.filenames, "r").read().strip().split("\n")
    genLen = len(genomeList)

    if args.filenames and genLen >=2:
        combs = combinations(genomeList,2)
        combLen = float(len(list(combinations(genomeList,2))))
        xcutoff = args.cre_cutoff
        ycutoff = args.acf_cutoff
        outdata_cre = []
        acf_res = []
        pool = mp.Pool(processes=args.nproc)

        if args.closely_related:
            kl = args.klarge
            acf_kmax2 = range(4, kl+1)
            pool1 = mp.Pool(processes=args.nproc)
            for l in range(genLen):
                outdata_cre.append(pool1.apply_async(kopt_cre, args=(l, kl, genomeList,args)))
            outdata_cre = list(chain(*[result.get() for result in outdata_cre]))
            cre = [x[1] for x in outdata_cre]
            if len(cre) != 0:
                maxval = xcutoff * float(max(cre))
            else:
                maxval = 0
            result = [x[0] for x in outdata_cre if x[1] < maxval]
            result2 = filter(result, genLen * 0)

            pool1.close()
            pool1.join()
            time.sleep(10)

            pool2 = mp.Pool(processes=args.nproc)
            outdata_ofc = []
            if len(result2) != 0:
                ks = min(result2)
            else:
                ks = 4
            for kmer in range(ks, kl+1):
                outdata_ofc.append(pool2.apply_async(kopt_ofc, args=(genomeList, kmer,args)))
            outdata_ofc = list(chain(*[result.get() for result in outdata_ofc]))
            ofc = [x[1] for x in outdata_ofc]
            if len(ofc) != 0:
                ofc_max = max(ofc)
            else:
                ofc_max = 1e7

            kmax = [x for x,y in outdata_ofc if y==ofc_max]
            ofc_possible_kmers = [x[0] for x in outdata_ofc if x[0] >= kmax[0]]
            pool2.close()
            pool2.join()

        else:
            for p in combs:
                acf_res.append(pool.apply_async(kopt_acf, args=(p,args)))
            acf_res = list(chain(*[result.get() for result in acf_res]))
            k_acf = [x[0] for x in acf_res]
            acf = [x[1] for x in acf_res]
            maxacf = ycutoff * float(max(acf))
            acf_kmax = [x[0] for x in acf_res if x[1] >= maxacf]
            acf_kmax2 = filter(acf_kmax, combLen * 0 )

            if len(acf_kmax2) == 0:
                kl = args.klarge
            else:
                kl = max(acf_kmax2)

            pool.close()
            pool.join()
            time.sleep(10)

            pool1 = mp.Pool(processes=args.nproc)
            for l in range(genLen):
                outdata_cre.append(pool1.apply_async(kopt_cre, args=(l, kl, genomeList,args)))
            outdata_cre = list(chain(*[result.get() for result in outdata_cre]))
            cre = [x[1] for x in outdata_cre]
            if len(cre) != 0:
                maxval = xcutoff * float(max(cre))
            else:
                maxval = 0
            result = [x[0] for x in outdata_cre if x[1] < maxval]
            result2 = filter(result, genLen * 0)

            pool1.close()
            pool1.join()
            time.sleep(10)
            
            pool2 = mp.Pool(processes=args.nproc)
            outdata_ofc = []
            if len(result2) != 0:
                ks = min(result2)
            else:
                ks = 4
            for kmer in range(ks, kl+1):
                outdata_ofc.append(pool2.apply_async(kopt_ofc, args=(genomeList, kmer,args)))
            outdata_ofc = list(chain(*[result.get() for result in outdata_ofc]))
            ofc = [x[1] for x in outdata_ofc]
            if len(ofc) != 0:
                ofc_max = max(ofc)
            else:
                ofc_max = 1e7

            kmax = [x for x,y in outdata_ofc if y==ofc_max]
            ofc_possible_kmers = [x[0] for x in outdata_ofc if x[0] >= kmax[0]]

            pool2.close()
            pool2.join()

    else:
        return print('list of filenames must be more than 1')

    opt = list(set(result2) & set(acf_kmax2) & set(ofc_possible_kmers))

    if args.output is not None and len(opt)>0:
        with open(args.output, 'w') as ofhandle:
            print("The recommended choice of kmer is: ", min(opt), file=ofhandle)
    
    elif args.output is not None and len(opt)==0:
        with open(args.output, 'w') as ofhandle:
            print(
                ("The recommended choice of kmer does not lie within the given interval "
                 "from 4 to the largest kmer length you selected"), 
                file=ofhandle
            )
    
    elif len(opt) == 0:
        print(("The recommended choice of kmer does not lie within the given interval "
               "from 4 to the largest kmer length you selected"))

    else:
        print("The recommended choice of kmer is: ", min(opt))


if __name__ == "__main__":
    run(sys.argv)
