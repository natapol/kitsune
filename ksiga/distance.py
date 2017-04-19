# -*- coding: utf-8 -*-

""" Collection of distance calculation function
 The distance matrix return as a distance array which constructs from column major lower triangle
 [0 1 2 3]
 [1 0 4 5] -> [0 1 2 3 4 5 6]
 [2 4 0 6]
 [3 5 6 0]
"""

import itertools
from itertools import zip_longest

import numpy as np
from sklearn.metrics import pairwise
from scipy.spatial.distance import squareform
from joblib import Parallel, delayed

import ksiga.mmath as mmath

def jensen_distance(m, n_thread=2):
    """ Calculate distance matrix using jensen distance

    Args:
        sm (matrix/sparse matrix): A normalize matrix.

    Returns: distance array 

    """
    rowNum = m.shape[0]

    if n_thread == 1:
        result = []
        for i, j in itertools.combinations(range(rowNum), r=2):
            iRow = m[i]
            jRow = m[j]
            distance = mmath.sparse_js_distance(iRow, jRow)
            result.append(distance)

        result = np.array(result)
    else: # Parallel version
        with Parallel(n_jobs=n_thread, backend="threading") as pool:
            generator = ((i, j) for i, j in itertools.combinations(range(rowNum), r=2))
            chunks = _grouper(100000, generator)
            result = pool(delayed(_calculate_chunk)(chunk, m)
                        for chunk in chunks)

    return np.squeeze(result)

#
# Helper functions
#
def _calculate_chunk(chunks, matrix):
    """ Calculate js distance in chunk.
    """

    results = []
    fchunks = filter(None, chunks)
    for i, j in fchunks:
        results.append(mmath.sparse_js_distance(matrix[i], matrix[j]))
    return results


def _grouper(n, iterable, padvalue=None):
    "_grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

#
# End of Helper functions
#


def cosine_distance(m, n_thread=1):
    """Calculate distance matrix using cosine

    Args:
        m (matrix/sparse matrix): TODO

    Returns: TODO

    """
    dm = pairwise.cosine_distances(m)
    da = squareform(dm, checks=False)
    return da


def euclidian_distance(m, n_thread=1):
    """Calculate distance matrix using euclidian distance

    Args:
        sm (matrix/sparse matrix): A normalized matrix.

    Returns: TODO

    """
    dm = pairwise.euclidean_distances(m)
    da = squareform(dm, checks=False)
    return da

def jaccard_distance(csr, n_thread=1, epsilon=1):
    """Calculate distance matrix using jaccard distance

    Args:
        sm (matrix/sparse matrix): A normalized matrix.

    Returns: TODO

    """
    import scipy.sparse as sp
    assert(0 < epsilon <= 1)

    csr = sp.csr_matrix(csr).astype(bool).astype(int)

    csr_rownnz = csr.getnnz(axis=1)
    intrsct = csr.dot(csr.T)

    nnz_i = np.repeat(csr_rownnz, intrsct.getnnz(axis=1))
    unions = nnz_i + csr_rownnz[intrsct.indices] - intrsct.data
    dists = 1.0 - intrsct.data / unions

    mask = (dists > 0) & (dists <= epsilon)
    data = dists[mask]
    indices = intrsct.indices[mask]

    rownnz = np.add.reduceat(mask, intrsct.indptr[:-1])
    indptr = np.r_[0, np.cumsum(rownnz)]

    dm = sp.csr_matrix((data, indices, indptr), intrsct.shape).toarray()
    da = squareform(dm, checks=False)
    return da


DISTANCE_FUNCTION = {
        "jensen" : jensen_distance,
        "cosine" : cosine_distance,
        "euclid" : euclidian_distance,
        "jaccard": jaccard_distance
        }
