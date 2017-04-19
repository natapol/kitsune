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
    return searchsorted(idx, index)


def searchsorted(haystack, needle):
    """ Search for a needle in haystack.
        np.searchsorted return an insert site. This one really do a searching.
        ** haystack have to be sorted.
    """
    idx = np.searchsorted(haystack, needle)
    mask = idx < haystack.size
    mask[mask] = haystack[idx[mask]] == needle[mask]
    idx = idx[mask]
    return idx
