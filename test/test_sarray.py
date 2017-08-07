# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sps

from ksiga import sparse_util as su 

def test_has_indices():
    # As condense as it can be
    data = [1,5,8,10,20,25]
    ind = [0,1,2,3,4,5]
    indptr = [0,5]

    m = sps.csr_matrix((data, ind, indptr), shape=(1,5))

    for i, idx in enumerate(ind):
        res = su.has_indices(m, idx)
        assert i == res

    # Some more test matrix
    data = [1,6,8,10,20,50,100]
    ind = [0,20,40,60,80,100, 120]
    indptr = [0,7]

    m = sps.csr_matrix((data, ind, indptr), shape=(1,500))

    for i, idx in enumerate(ind):
        res = su.has_indices(m, idx)
        assert i == res


def test_sort_indice():
    """ Test if an indice of sprase array still sort after addition.
        The reason of this test because our method rely on the sorted indice.
    """

    # I am lazy.
    for i in range(50):
        a = sps.rand(1,50000, 0.2).tocsr()
        b = sps.rand(1,50000, 0.2).tocsr()
        assert a.has_sorted_indices == 1
        assert b.has_sorted_indices == 1
        assert (a + b).has_sorted_indices == 1
