# -*- coding: utf-8 -*-

""" Store and return k-mer
"""

import os
import re
import datetime

import Bio.SeqIO as SeqIO
import numpy as np
import h5py
import scipy.sparse as sps

import ksiga
from ksiga.kmerutil import kmer_hash_emit_rolling


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


def split_sequence(seqioH, ksize):
    """ Split and clean sequence into piece to work with.
    """
    notAllowReg = re.compile("[^ATCG]")
    for seq in seqioH:
        rawStr = str(seq.seq).upper()
        rawSeqs = notAllowReg.split(rawStr)
        for rawSeq in rawSeqs:
            if len(rawSeq) < ksize:
                continue
            yield rawSeq

def build_csr_matrix_from_fasta(fh, ksize):
    # Check total number of base
    handleCount = SeqIO.parse(fh, "fasta")
    count = 1
    # Init array first.
    for seq in handleCount:
        count += len(seq)
    handleCount.close()
    emits = np.zeros(count, dtype=np.int64)  # Initialize array to hold data.

    # Reinit iterator again, now to count
    fh.seek(0, 0)
    handle = SeqIO.parse(fh, "fasta")
    locpointer = 0
    for seq in split_sequence(handle, ksize):
        for hsh in kmer_hash_emit_rolling(seq, ksize):
            emits[locpointer] = hsh
            locpointer += 1

    indices, data = np.unique(emits[0:locpointer], return_counts=True)

    return indices, data

def build_signature(fasta, ksize, store, force):
    """ Build signature file from fasta.

    Args:
        fasta (fh): filehandle for fasta file
        ksize (int): kmer size
        store (str): filename
        force (bool): Rewrite

    Returns: Write a signature into file

    """

    if force:
        store = h5py.File(store, "w")
    else:
        store = h5py.File(store, "a")
    sigName = "/kmer/{size}".format(size=ksize)
    if sigName in store:
        raise KeyError("Kmer is already create")

    group = store.create_group(sigName)
    indice, data = build_csr_matrix_from_fasta(fasta, ksize)
    colSize = 4 ** (ksize)
    group.create_dataset("indices", compression="gzip", compression_opts=5, data=indice, dtype='int64')
    group.create_dataset("data", compression="gzip", compression_opts=5, data=data, dtype='int64')
    group.create_dataset("shape", data=(1, colSize))
    # Create metadata
    group.attrs["filename"] = fasta.name
    group.attrs["run_date"] = datetime.datetime.now().isoformat()
    group.attrs["run_version"] = ksiga.__version__
    group.attrs["run_status"] = "finished"