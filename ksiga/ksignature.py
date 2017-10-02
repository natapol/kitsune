# -*- coding: utf-8 -*-

""" Store and return k-mer
"""

import os

import numpy as np
import h5py
import scipy.sparse as sps


class KmerSignature(object):

    """Store and represent signature of kmer. It is a wrapper around HDF5 format, where kmer is store"""

    # How to store sparse matrix.
    KMER_TP = "/kmer/{size}"
    KMERD_TP = "/kmer/{size}/data"
    KMERI_TP = "/kmer/{size}/indices"
    KMERS_TP = "/kmer/{size}/shape"

    def __init__(self, pointer):
        """TODO: Docstring for function.

        Args:
            pointer (str): file name to store index in.

        Returns: TODO

        """

        if not os.path.exists(pointer):
            raise FileNotFoundError("Cannot find an index file")

        self._pointer = h5py.File(pointer, "r")

    def get_kmer(self, ksize):
        """ Get kmer from pointer.

        Args:
            ksize (int): size of kmer to calculate.

        Returns: Sparse matrix of one row

        # TODO Make a optional to return what each column kmer is?

        """

        if not self.has_kmer(ksize):
            raise IndexError("K-mer size = {}, have not yet index.".format(ksize))

        indices = self._pointer[self.KMERI_TP.format(size=ksize)]
        data = self._pointer[self.KMERD_TP.format(size=ksize)]
        shape = self._pointer[self.KMERS_TP.format(size=ksize)]
        # It is trivial, but check this just in case.
        assert indices.shape == data.shape
        return sps.csr_matrix((data, indices, [0, data.shape[0]]), shape=shape)

    def has_kmer(self, ksize):
        """TODO: Check if k-mer of length is indexed.

        Args:
            ksize (int): TODO

        Returns: TODO

        """

        access = self.KMER_TP.format(size=ksize)
        return access in self._pointer

    def list_all_kmer(self):
        """ List all kmer
        Returns: list of kmer

        """
        return list(self._pointer["kmer"])

    @classmethod
    def fromfilename(cls, filename):
        """ Initialize data from filename"""
        return cls(filename)
        

    @classmethod
    def fromstore(cls, store):
        """ Initialize data from filename"""
        return store


def rebuild_sparse_matrix(stores, ksize):
    """ Rebuild sparse matrix from list of stores. The implementation
        is obvious (if not, read csr_matrix doc on scipy).

    Args:
        stores (TODO): TODO
        ksize (int): Size of kmer to calculate

    Returns: TODO

    """
    # Initialize list for building a sparse matrix.
    data = []
    indices = []
    indptr = []

    indptr.append(0)

    for store in stores:
        array = KmerSignature(store).get_kmer(ksize)
        data.append(array.data)
        indices.append(array.indices)
        indptr.append(indptr[-1] + array.data.shape[0])

    data = np.concatenate(data).astype(np.int64)
    indices = np.concatenate(indices).astype(np.int64)

    rowNum = len(stores)
    colNum = 4 ** ksize
    shape = (rowNum, colNum)

    return sps.csr_matrix((data, indices, indptr), shape=shape)
