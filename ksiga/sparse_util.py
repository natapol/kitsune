# -*- coding: utf-8 -*-

""" 
"""

import numpy as np

# TODO: Refactor operation into seperate class call SparseArray.
# TODO: This one is pretty much the same with sortedsearch. REFACTOR it.
def has_indices(sparseArray, index):
    """ Search if sparse array has index.
        Return: Position of index if found, empty array if not.
    """
    # Search in sparse array
    ind = sparseArray.indices
    ret = searchsorted(ind, np.array([index]))
    if list(ret):  #  Numpy's behavior is a bit different to list.
        return ret[0]
    else:
        return ind.shape[0]


# def has_indices(sparseArray, index):
    # """ Search if sparse array has index.
        # Return: Position of index if found, size of an array if not found.
    # """
    # # Search in sparse array
    # ind = sparseArray.indices

    # # Must be sorted
    # insert = np.searchsorted(ind, index)
    # try:
        # if ind[insert] == index:
            # return insert
        # else:
            # return ind.shape[0]
    # except IndexError:
        # return ind.shape[0]

def searchsorted(haystack, needle):
    """ Search for a needle in haystack.
        np.searchsorted return an insert site. This one really do a searching.
        ** haystack have to be sorted.
    Args:
        haystack (np.array): TODO
        needle (np.array): TODO

    Returns: TODO
        
    """
    idx = np.searchsorted(haystack, needle)
    mask = idx < haystack.size
    mask[mask] = haystack[idx[mask]] == needle[mask]
    idx = idx[mask]
    return idx
