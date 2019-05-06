"""
.. module:: cre
   :platform: Unix, MacOSX
   :synopsis: module for Cumulative Relative Entropy (CRE) calculation

.. moduleauthor:: Natapol Pornputtapong <natapol.p@chula.ac.th>


"""
from . import kitsunejf as jf
import math
from tqdm import tqdm

def cal_cre(fsa, kend, kfrom=4, **karg):
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
        left = a1[key[0:-1]]
        right = a1[key[1:]]
        below = a2[key[1:-1]]
        if 0 not in (left, right, below):
            expectf = left * right / below
            result += max(0, realf / rfactor * (math.log(realf / expectf, 2) + lfactor))
    return result

if __name__ == "__main__":
    print(cal_cre('./examples/ASM170810v1_genomic.fna', 25))