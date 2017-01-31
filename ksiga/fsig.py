# -*- coding: utf-8 -*-

""" Store and calculate k-mer.
    Also calculate various statistic of kmer (1).
    (1) Viral Phylogenomics using an alignment-free method: A three step approach to dtermine optimal length of k-mer
"""

import os

import scipy.sparse as sps
import h5py
import numpy as np

from ksiga import kmer
from ksiga import logutil


class KmerSignature(object):

    """Store and represent signature of kmer. Currently is a wrapper around HDF5 format"""

    # Folder, Data, Indices and Shape
    KMER_TP = "/kmer/{size}"
    KMERD_TP = "/kmer/{size}/data"
    KMERI_TP = "/kmer/{size}/indices"
    KMERS_TP = "/kmer/{size}/shape"

    def __init__(self, pointer):
        """TODO: Docstring for function.

        Args:
            pointer (folder): TODO

        Returns: TODO

        """

        if not os.path.exists(pointer):
            raise FileNotFoundError("Sumting wong")

        self._pointer = h5py.File(pointer, "r")


    def get_kmer(self, ksize):
        """ Get kmer from pointer.

        Args:
            ksize (int): size of kmer to calculate.

        Returns: Sparse matrix of one row

        # TODO Make a optional to return what each column kmer is?

        """

        if not self.has_kmer(ksize):
            raise IndexError("No kmer")

        indices = self._pointer[self.KMERI_TP.format(size=ksize)]
        data = self._pointer[self.KMERD_TP.format(size=ksize)]
        shape = self._pointer[self.KMERS_TP.format(size=ksize)]

        # Obvious, but just in case.
        assert indices.shape == data.shape
        return sps.csr_matrix((data, indices, [0, data.shape[0]]), shape=shape)


    def has_kmer(self, ksize):
        """TODO: Docstring for has_kmer.

        Args:
            ksize (TODO): TODO

        Returns: TODO

        """

        access = self.KMER_TP.format(size=ksize)
        return access in self._pointer


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
def has_indices(sparseArray, index):
    """ Search if sparse array has index.
        Return size of array if not found ()
    """
    # Search in sparse array
    ind = sparseArray.indices

    # If sorted
    insert = np.searchsorted(ind, index)
    try:
        if ind[insert] == index:
            return insert
        else:
            return ind.shape[0]
    except IndexError:
        return ind.shape[0]

def calculate_relative_entropy(store, ksize):
    """ Calculate the relative entropy (obs vs expected)
        The equation is (Something in latex here)
        Expected (ABCD) =   (obs(ABC) * obs(BCD)) / obs(BC)
    """
    # Check if exists
    store = KmerSignature(store)

    # Sparse array of ONE row, so it has shape = (1, col)
    array0 = store.get_kmer(ksize)
    array1 = store.get_kmer(ksize-1)
    array2 = store.get_kmer(ksize-2)

    merGen = kmer.generateMers(ksize)

    genLoc0 = kmer.create_kmer_loc_fn(ksize)
    genLoc1 = kmer.create_kmer_loc_fn(ksize-1)
    genLoc2 = kmer.create_kmer_loc_fn(ksize-2)

    ARes = []
    FRes = []
    BRes = []
    MRes = []

    # Calculate final normalization factor. Delegate the normalization to the last step (Arithmatric).

    # TODO: May be it would be easier to do it in dense array? Dense array take a lot of space
    # but it does not need to do all of these magic.
    # TODO: Instead of by mer, just iterate the biggest kmer array.
    for mer in merGen:
        merFront = mer[0:-1]
        merBack = mer[1:]
        merMiddle = mer[1:-1]

        locL = genLoc0(mer)
        # Normal csr procedure is to convert to coo first but we can skip all of that since
        # 1. We make sure it sorted 2. We don't need Y-axis location.
        idxA = has_indices(array0, locL)

        if idxA == array0.indices.shape[0]:
            # This kmer does not have any obs. skip it.
            continue

        # logutil.notify(mer)
        # logutil.notify(locL)
        # logutil.notify(idxA)
        # logutil.notify(array0.indices[0:10])
        locF = genLoc1(merFront)
        locB = genLoc1(merBack)
        locM = genLoc2(merMiddle)

        # Find location of left mer, right mer, and middle mer
        idxF = has_indices(array1, locF)
        idxB = has_indices(array1, locB)
        idxM = has_indices(array2, locM)

        # For debugging.
        if (idxF == array1.indices.shape[0]):
            raise IndexError("Left not found")
        if (idxB == array1.indices.shape[0]):
            raise IndexError("Right not found")
        if (idxM == array2.indices.shape[0]):
            raise IndexError("Middle not found")

        # All, Front, Back, Middle
        countA = array0.data[idxA]
        countF = array1.data[idxF]
        countB = array1.data[idxB]
        countM = array2.data[idxM]

        ARes.append(countA)
        FRes.append(countF)
        BRes.append(countB)
        MRes.append(countM)

    # Calculate.
    ARes = np.array(ARes)  # Obs
    FRes = np.array(FRes)  # Front
    BRes = np.array(BRes)  # Back
    MRes = np.array(MRes)  # Middle

    # Factor version
    norm0 = array0.data.sum()
    norm1 = array1.data.sum()
    norm2 = array2.data.sum()
    expectation = (FRes * BRes) / MRes
    observation = ARes
    normFactor = (norm1 ** 2) / (norm2 * norm0)
    rhs = np.log2(observation / expectation) + np.log2(normFactor)
    lhs = ARes / norm0
    relativeEntropy = (lhs * rhs).sum()
     
    return relativeEntropy

def calculate_common_feature(store1, store2, ksize):
    """ Calculate the common feature between pair of sequences

    Args:
        store1 (list(str)): Storage
        store2 (list(str)): Storage
        ksize (int): k-mer size.

    Returns: TODO

    """
    # Build sparse matrix.
    store1 = h5py.File(store1, "r")
    store2 = h5py.File(store2, "r")

    sigName = "/kmer/{size}".format(size=ksize)

    if (sigName not in store1) or (sigName not in store2):
        raise IndexError("Not kmer")

def calculate_average_common_feature(csrMatrix, ksize):
    """ Calculate an average common features from sparse matrix.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    N = csrMatrix.shape[0]


def rebuild_sparse_matrix(stores, ksize):
    """TODO: Docstring for rebuild_sparse_matrix.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """

    # data 
    # indices
    # indptr
    # for store in stores:
        


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
    shapeSize = 4 ** (ksize)
    group.create_dataset("indices", compression="gzip", compression_opts=5, data=indice)
    group.create_dataset("data", compression="gzip", compression_opts=5, data=data)
    group.create_dataset("shape", data=(1, shapeSize))
