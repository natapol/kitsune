# -*- coding: utf-8 -*-

""" Store and calculate k-mer.
    Also calculate various statistic of kmer (1).
    (1) Viral Phylogenomics using an alignment-free method: A three step approach to dtermine optimal length of k-mer
"""

import scipy.sparse as sps
import h5py
import numpy as np

from ksiga.ksignature import KmerSignature
from ksiga import kmer
from ksiga import logutil

# TODO: Refactor operation into seperate class call SparseArray.
# TODO: This one is pretty much the same with sortedsearch. REFACTOR it.
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


def sortedsearch(npArray, vals):
    """ Search for sorted array
        Unlike numpy's implementation, this return an vals that found

    Args:
        npArray (TODO): TODO
        vals (TODO): TODO

    Returns: numpy array of vals that found in npArray

    """
    idx = np.searchsorted(npArray, vals)
    mask_idx = idx == npArray.shape[0]
    idx[mask_idx] = npArray.shape[0] - 1
    valsIndex = npArray[idx] == vals
    return vals[valsIndex]


def calculate_relative_entropy(store, ksize):
    """ Calculate the relative entropy (obs vs expected)
        The equation is (Something in latex here)
        Expected (ABCD) =  (obs(ABC) * obs(BCD)) / obs(BC)
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


def calculate_average_common_feature(stores, ksize):
    """ Calculate an average common features from sparse matrix.

    Args:
        stores (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    import itertools

    csr_matrix = rebuild_sparse_matrix(stores, ksize)
    rowNum = csr_matrix.shape[0]

    counts = []
    norm = csr_matrix.shape[0] - 1

    for i, j in itertools.combinations(range(rowNum), r=2):
        iRow = csr_matrix[i]
        jRow = csr_matrix[j]
        # Find column (same kmer) that contain in both row. Good thing that we already sorted it.
        found = sortedsearch(iRow.indices, jRow.indices)
        counts.append(found.shape[0])

    return sum(counts) / norm


def calculate_obsff(store, ksize):
    """ Calculate an average common features from sparse matrix.

    Args:
        store (TODO): TODO
        ksize (TODO): TODO

    Returns: TODO

    """
    store = KmerSignature(store)
    # Sparse array of ONE row, so it has shape = (1, col)
    array = store.get_kmer(ksize)
    count = array.indices.shape[0]

    expect = 4 ** ksize
    return count / ksize


def rebuild_sparse_matrix(stores, ksize):
    """ Rebuild sparse matrix from list of stores

    Args:
        stores (TODO): TODO
        ksize (int): Size of kmer to calculate

    Returns: TODO

    """

    # I don't know why scipy loves to convert sparse matrix to COO first....
    data = []
    indices = []
    indptr = []

    indptr.append(0)

    for store in stores:
        array = KmerSignature(store).get_kmer(ksize)
        data.append(array.data)
        indices.append(array.indices)
        indptr.append(indptr[-1] + array.data.shape[0])

    data = np.concatenate(data)
    indices = np.concatenate(indices)

    rowNum = len(stores)
    colNum = 4 ** ksize
    shape = (rowNum, colNum)

    return sps.csr_matrix((data, indices, indptr), shape=shape)


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
