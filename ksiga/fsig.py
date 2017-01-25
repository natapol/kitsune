# -*- coding: utf-8 -*-

""" Store and release kmer freq.
"""

import h5py
import numpy as np

from ksiga import kmer


class KmerSignature(object):

    """Store and represent signature of kmer. Currently is a wrapper around HDF5 format"""

    def __init__(self, pointer):
        """TODO: Docstring for function.

        Args:
            pointer (folder): TODO

        Returns: TODO

        """
        self._pointer = pointer

        if not os.path.exists(self._pointer):
            raise FileNotFoundError("Sumting wong")

    def get_kmer(self, ksize):
        """ Get kmer from pointer.

        Args:
            ksize (int): 

        Returns: TODO

        """

        pass

    def has_kmer(self, ksize):
        """TODO: Docstring for has_kmer.

        Args:
            ksize (TODO): TODO

        Returns: TODO

        """
        pass

    def list_all(self):
        """ List all index here
        Returns: TODO

        """
        pass

    # # @lru_cache
    # def kmer_location(kmer):
        # """ Calculate kmer location to store in array
        # """
        # encMap = {"A":1, "T":2, "C":3, "G":4}

        # code = 0
        # for ch in kmer:
            # code *= 4
            # code += encMap[ch]

        # return code

# TODO: Refactor operation into seperate class call SparseArray.

def calculate_relative_entropy(store, ksize):
    """ Calculate the relative entropy (obs vs expected)
        The equation is (Something in latex here)
    """
    # Check if exists
    sigL = "/kmer/{size}".format(size=ksize)
    sigS = "/kmer/{size}".format(size=ksize-1)

    if (sigL not in store) or (sigS not in store):
        raise IndexError("Specific kmer lenght haven't been calculated yet")

    # Sparse array of ONE row
    arrayL = store.get_kmer(sigL)
    arrayS = store.get_kmer(sigS)

    merGen = kmer.generateMers(ksize)

    genLocL = kmer.create_kmer_loc_fn(ksize)
    genLocS = kmer.create_kmer_loc_fn(ksize-1)

    for mer in merGen:
        merFront = mer[0:-1]
        merBack = mer[1:]
        locL = genLocL(mer)
        # Normal csr procedure is to convert to coo first but we can skip all of that since
        # 1. We make sure it sorted 2. We don't need Y-axis location.
        idx = np.searchsorted(arrayL, locL.indices)
        if idx == arrayL.shape[1]:
            # This kmer does not have any obs. skip it.
            continue

        locF = genLocS(merFront)
        locB = genLocS(merBack)
        idxF = np.searchsorted(merFront, arrayS.indices)
        idxB = np.searchsorted(merBack, arrayS.indices)

        if (idxF == arrayS.shape[1]) or (idxB == arrayS.shape[1]):
            raise IndexError("We are a big trouble.")

        arrayS.data[idxF]
        arrayS.data[idxB]

        # arrayL.indices
        # locFront = genLocSmall()
        # locBack = gen

def calculate_common_feature(store1, store2, ksize):
    """ Calculate the common feature between pair of sequences

    Args:
        store1 (str): TODO
        store2 (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    store1 = h5py.File(store1, 'r')
    store2 = h5py.File(store2, 'r')

    sigName = "/kmer/{size}".format(size=ksize)

    if (sigName not in store1) or (sigName not in store2):
        raise IndexError("Not kmer")

def build_signature(fasta, ksize, store):
    """ Build signature file from fasta.

    Args:
        fasta (str): filename
        ksize (int): kmer size
        store (str): filename

    Returns: Write a signature into file

    """

    store = h5py.File(store, "a")
    sigName = "/kmer/{size}".format(size=ksize)
    if sigName in store:
        raise KeyError("Already existed")

    group = store.create_group(sigName)
    indice, data = kmer.build_csr_matrix_from_fasta(fasta, ksize)
    group.create_dataset("indices", compression="gzip", compression_opts=5, data=indice)
    group.create_dataset("data", compression="gzip", compression_opts=5, data=data)
