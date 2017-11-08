# -*- coding: utf-8 -*-

import pytest

import numpy as np
import scipy.sparse as sp
import ksiga.sparse_util as su


def test_searchsorted_all2all():
    """ Shouldn't miss any
    """
    for i in range(20):
        sample = np.sort(np.random.choice(np.arange(0,1000), 100, replace=False))
        expected = np.arange(0,100)
        answer = su.searchsorted(sample, sample)
        print(answer)
        assert np.array_equal(answer, expected)

def test_searchsorted_subset():
    """ When query is a subset of the target
    """
    for i in range(20):
        sample = np.sort(np.random.choice(np.arange(0,1000), 100, replace=False))
        remove = np.random.choice(np.arange(0,100), 10, replace=False)
        subsample = np.delete(sample, remove)
        expected = np.delete(np.arange(0, 100), remove)
        answer = su.searchsorted(sample, subsample)
        assert np.array_equal(answer, expected)

def _random_generate_csr(size):
    pass

def test_sparse_array_match():
    csr_array1 = sp.random(1, 50, density=0.9, format="csr")
    csr_array2 = sp.random(1, 50, density=0.9, format="csr")

    # Check with intersect.
    answer = set(csr_array1.indices).intersection(csr_array2.indices)

    

def test_searchsorted():
    test_1 = np.array([1,3,5,7,9])
    test_2 = np.array([3,6,9])
    assert np.array_equal(su.searchsorted(test_1, test_2), [1,4])
    assert np.array_equal(su.searchsorted(test_2, test_1), [0,2])