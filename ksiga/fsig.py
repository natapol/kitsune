# -*- coding: utf-8 -*-

""" Store and calculate various statistic of kmer (1).
    (1) Viral Phylogenomics using an alignment-free method: A three step approach to determine optimal length of k-mer
"""

import scipy.sparse as sps
from scipy.interpolate import interp1d
import h5py
import numpy as np

import ksiga.ksignature as ksignature
from ksiga.ksignature import KmerSignature
import ksiga.sparse_util as su
from ksiga import kmer
from ksiga import logutil


def calculate_relative_entropy(indexFilename, ksize):
    """ Calculate the relative entropy (obs vs expected)
        The equation is (Something in latex here)
        Expected (ABCD) =  (obs(ABC) * obs(BCD)) / obs(BC)
    """
    # Check if exists
    store = KmerSignature(indexFilename)
    # Sparse array of ONE row, so it has shape = (1, col)
    array0 = store.get_kmer(ksize)
    array1 = store.get_kmer(ksize-1)
    array2 = store.get_kmer(ksize-2)

    genLoc1 = kmer.create_kmer_loc_fn(ksize-1)
    genLoc2 = kmer.create_kmer_loc_fn(ksize-2)

    ARes = []
    FRes = []
    BRes = []
    MRes = []
    # Calculate final normalization factor. Delegate the normalization to the last step (Arithmatric).
    # TODO: Collect everything and calculate in vectorize manner.?
    for (idxA, locL) in enumerate(array0.indices):
        mer = kmer.decode(locL, ksize)
        merFront = mer[0:-1]
        merBack = mer[1:]
        merMiddle = mer[1:-1]
        # Find location of left mer, right mer, and middle mer
        locF = genLoc1(merFront)
        locB = genLoc1(merBack)
        locM = genLoc2(merMiddle)
        # TODO: There should be an easy and very effieicient way to map ATTT -> ATT, TTT
        idxF = su.has_indices(array1, locF)
        idxB = su.has_indices(array1, locB)
        idxM = su.has_indices(array2, locM)
        # For debugging. This shouldn't be happened
        if __debug__:
            if idxF == array1.indices.shape[0]:
                raise IndexError("Left not found")
            if idxB == array1.indices.shape[0]:
                raise IndexError("Right not found")
            if idxM == array2.indices.shape[0]:
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

    ## Version which follow a formular more closely.
    ## Roughly the same speed with the reduce version above.
    # norm0 = array0.data.sum()
    # norm1 = array1.data.sum()
    # norm2 = array2.data.sum()
    # expectation = (FRes/norm1) * (BRes/norm1) / (MRes/norm2)
    # observation = ARes / norm0
    # relativeEntropy = (observation * np.log2(observation/expectation)).sum()

    return relativeEntropy


def calculate_cre(indexFilename, start_k, stop_k):
    """ Calculate CRE and 

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """

    #  Check that all k-mer is already indexed.
    store = KmerSignature(indexFilename)
    for k in range(start_k - 2, stop_k + 1):
        pass
        # raise IndexError()

    relEntropy = []
    for k in range(start_k+1, stop_k+1):  # You don't need to calculate the first one (because it is equal to sum).
        relEntropy.append(calculate_relative_entropy(indexFilename, k))

    relEntropy = np.array(relEntropy)
    #  This shouldn't have much impact, but just in case where number goes below 0.
    relEntropy = np.clip(relEntropy, 0, np.inf)
    maximumVal = relEntropy.sum()
    cre = maximumVal - relEntropy.cumsum()
    cre = np.insert(cre, 0, maximumVal)
    kmerRange = np.arange(start_k, stop_k+1)
    suggestKmer = _find_yintercept(kmerRange, cre, 10)

    return (cre, suggestKmer)

def _find_yintercept(x, y, percent):
    """ Find an intercept of CRE graph. ( Log/Ln graph)?
    """

    cutoff = y.max() / percent
    #  Check if it is still fall in range.
    if y.min() > cutoff:
        raise NotImplementedError("Calculate the k-mer beyond index is not implemented yet.")
    #  Interpolate to get more resolution on k-mer.
    fi = interp1d(x, y, kind="cubic")
    xn = np.linspace(x[0], x[-1], 500)
    yn = fi(xn)
    #  Check for an intercept
    idx = np.argwhere(np.diff(np.sign(yn - cutoff)) != 0).reshape(-1)
    xIntercept = xn[idx[0]]
    kmer = int(round(xIntercept)) #  Kmer
    return kmer


def calculate_average_common_feature(stores, ksize):
    """ Calculate an average common features from sparse matrix.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """

    csr_matrix = rebuild_sparse_matrix(stores, ksize)
    rowNum = csr_matrix.shape[0]

    vals = []
    norm = csr_matrix.shape[0] - 1

    # TODO: Try C * tranpose(C) when all data are convert to 1.
    # That should somewhat faster?
    # 1. http://stackoverflow.com/questions/24566633/which-is-the-best-way-to-multiply-a-large-and-sparse-matrix-with-its-transpose
    # 2. C.dot(C.transpose)
    csr_matrix.data = np.ones(csr_matrix.data.shape[0], np.int64)
    for i in range(rowNum):
        val = 0
        # val = csr_matrix.dot(csr_matrix[i].transpose()).sum()  # Need to work on this later
        for j in range(rowNum):
            if i == j:
                continue
            iRow = csr_matrix[i]
            jRow = csr_matrix[j]

            found = su.searchsorted(iRow.indices, jRow.indices)
            val += found.shape[0]

        vals.append(val)

    result = np.array(vals) / norm

    return result


def calculate_acf(stores, start_k, stop_k):
    
    results = []

    for ksize in range(start_k, stop_k + 1):
        val = calculate_average_common_feature(stores, ksize)
        val = val[:, np.newaxis]
        # Turn into 
        results.append(val)

    results = np.hstack(results)
    # For each row, calculate the k-mer
    for result in results:
        kmer = _find_yintercept(np.arange(start_k, stop_k+1), result, 10)
        print(kmer)


def calculate_obsff(stores, ksize):
    """ Calculate an observe and expect.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    csr_matrix = rebuild_sparse_matrix(stores, ksize)
    # Probalbility that kmers exist.
    norm = csr_matrix.sum()
    prob = np.asarray(csr_matrix.sum(axis=0)).squeeze() / norm
    # Remove zero
    prob = prob[np.nonzero(prob)]
    # How many genome they occur
    csr_matrix.data = np.ones_like(csr_matrix.data)
    occurence = np.asarray(csr_matrix.sum(axis=0)).squeeze()
    occurence = occurence[np.nonzero(occurence)]  #  Remove zero
    sites = np.unique(csr_matrix.indices)  # Kmer string
    fn = np.vectorize(kmer.decode)
    kmerStr = fn(sites, ksize)

    return (prob, occurence, kmerStr)


def calculate_uniq_mer(stores, ksize):
    """ Calculate an observe and expect.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    csr_matrix = rebuild_sparse_matrix(stores, ksize)
    unique, count = np.unique(csr_matrix.indices, return_counts=True)

    allPossible = unique.shape[0]
    numberOfUnique = np.where(count == 1)[0].shape[0]

    return (allPossible, numberOfUnique)


rebuild_sparse_matrix = ksignature.rebuild_sparse_matrix

# def rebuild_sparse_matrix(stores, ksize):
    # """ Rebuild sparse matrix from list of stores. The implementation
        # is obvious (if not, read csr_matrix doc on scipy).

    # Args:
        # stores (TODO): TODO
        # ksize (int): Size of kmer to calculate

    # Returns: TODO

    # """
    # # Initialize list for building a sparse matrix.
    # data = []
    # indices = []
    # indptr = []

    # indptr.append(0)

    # for store in stores:
        # array = KmerSignature(store).get_kmer(ksize)
        # data.append(array.data)
        # indices.append(array.indices)
        # indptr.append(indptr[-1] + array.data.shape[0])

    # data = np.concatenate(data)
    # indices = np.concatenate(indices)

    # rowNum = len(stores)
    # colNum = 4 ** ksize
    # shape = (rowNum, colNum)

    # return sps.csr_matrix((data, indices, indptr), shape=shape)


def build_signature(fasta, ksize, store):
    """ Build signature file from fasta.

    Args:
        fasta (fh): filehandle for fasta file
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
    colSize = 4 ** (ksize)
    group.create_dataset("indices", compression="gzip", compression_opts=5, data=indice)
    group.create_dataset("data", compression="gzip", compression_opts=5, data=data)
    group.create_dataset("shape", data=(1, colSize))
