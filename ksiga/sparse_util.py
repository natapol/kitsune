# -*- coding: utf-8 -*-

""" Utility for sparse matrix
"""

import numpy as np

# TODO: Refactor operation into seperate class call SparseArray.
# TODO: This one is pretty much the same with sortedsearch. REFACTOR it.
def has_indices(sparseArray, index):
    """ Search if sparse array has index.
        Return: Position of index if found. Out of bound index if not found.
    """
    # Search in sparse array
    ind = sparseArray.indices
    ret = searchsorted(ind, np.array([index]))
    if list(ret):  #  Numpy's behavior is a bit different to list.
        return ret[0]
    else:
        return ind.shape[0]


def searchsorted(haystack, needle):
    """ Search for a needle in haystack.
        np.searchsorted return an insert site. This one really do a searching.
        ** haystack have to be sorted.
    Args:
        haystack (np.array): array that we want to search in. Must be sort and unique
        needle (np.array): query that want to search. Also sort and unique

    Returns: Array of location that found. Note that this SKIP the value that does not found.
        
    """
    idx = np.searchsorted(haystack, needle)
    mask = idx < haystack.size
    mask[mask] = haystack[idx[mask]] == needle[mask]
    idx = idx[mask]
    return idx

def get_matching(v1, v2):
    """ Get array of value where value is found in both vectors
    """
    pass