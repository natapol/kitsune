# -*- coding: utf-8 -*-

import numpy as np

def searchsorted(haystack, needle):
    """ ** BOTH ARRAYS ARE ASSUMED TO BE SORTED **
    """
    idx = np.searchsorted(haystack, needle)
    mask = idx < haystack.size
    mask[mask] = haystack[idx[mask]] == needle[mask]
    idx = idx[mask]
    return idx
