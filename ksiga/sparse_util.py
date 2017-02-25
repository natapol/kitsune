# -*- coding: utf-8 -*-

import numpy as np

# TODO: Refactor operation into seperate class call SparseArray.
# TODO: This one is pretty much the same with sortedsearch. REFACTOR it.
def has_indices(sparseArray, index):
    """ Search if sparse array has index.
        Return: Position of index if found, size of an array if not found.
    """
    # Search in sparse array
    ind = sparseArray.indices

    # Must be sorted
    insert = np.searchsorted(ind, index)
    try:
        if ind[insert] == index:
            return insert
        else:
            return ind.shape[0]
    except IndexError:
        return ind.shape[0]

def searchsorted(haystack, needle):
    """ ** BOTH ARRAYS ARE ASSUMED TO BE SORTED **
    """
    idx = np.searchsorted(haystack, needle)
    mask = idx < haystack.size
    mask[mask] = haystack[idx[mask]] == needle[mask]
    idx = idx[mask]
    return idx
