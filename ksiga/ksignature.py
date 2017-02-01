# -*- coding: utf-8 -*-

import os

import h5py
import scipy.sparse as sps

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
