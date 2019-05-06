"""
.. module:: acf
   :platform: Unix, MacOSX
   :synopsis: module for Average number of common features (ACF) calculation

.. moduleauthor:: Natapol Pornputtapong <natapol.p@chula.ac.th>


"""
from . import kitsunejf as jf
import math
from tqdm import tqdm

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
        result = {}
        for kmer in tqdm(kmers):
            keys_array = []
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

if __name__ == "__main__":
    print(cal_acf(['./examples/ASM170810v1_genomic.fna', './examples/ASM170810v1_genomic.fna'], [10,11]))
