#!/usr/bin/env python
from . import kitsunejf as jf
import math
from tqdm import tqdm
import numpy as np
import collections

def cal_cre(fsa, ksmall, klarge, **karg):
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

def cal_acf(fsas, klarge, ksmall, **karg):
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

    #for kmer in tqdm(range(klarge,ksmall-1, -1)):
    #    if kmer % 2 == 1:
    keys = []
    if kmer % 2 == 1:
        for fsa in fsas:
            keys.extend(list(jf.Kmercount(fsa, kmer, **karg).keys()))

        count_feature = list(collections.Counter((collections.Counter(keys).values())).values())
        lnp = np.log2(count_feature) - (2 * kmer)
        result[kmer] = np.sum(np.exp2(lnp) * lnp) * -1
    return result


if __name__ == "__main__":
    print(cal_cre('./examples/ASM170810v1_genomic.fna', 25))

