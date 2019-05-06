"""
.. module:: ofc
   :platform: Unix, MacOSX
   :synopsis: module for observed feature occurrences (OFC) calculation

.. moduleauthor:: Natapol Pornputtapong <natapol.p@chula.ac.th>


"""
from . import kitsunejf as jf
import math
import collections
import numpy as np
from tqdm import tqdm

def cal_ofc(fsas, kmers, **karg):
    """ Calculate shannon entropy of observed feature occurrences (OFC)
        ofc(l) = -sum(p ln p)
    Args:
        fsa a genome file
        kmers a list of kmer to calculate

    Returns: float shannon entropy

    """
    result = {}
    
    for kmer in tqdm(kmers):
        keys = []
        for fsa in fsas:
            keys.extend(list(jf.Kmercount(fsa, kmer, **karg).keys()))
        
        count_feature = list(collections.Counter((collections.Counter(keys).values())).values())
        lnp = np.log2(count_feature) - (kmer * 2)
        result[kmer] = np.sum(np.exp2(lnp) * lnp) * -1

    return result

if __name__ == "__main__":
    print(cal_ofc(['./examples/ASM170810v1_genomic.fna', './examples/ASM170810v1_genomic.fna'], [10,11,12]))
