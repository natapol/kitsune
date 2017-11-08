# -*- coding: utf-8 -*-

import pytest

import numpy as np
import ksiga.kmerutil as kmer


def _generate_sequence(l=500):
    # Generate a random genomic sequence
    # Loop is better, this is an excercise for
    # Broadcasting.
    p = np.array([0.3, 0.3, 0.2, 0.2])
    c = p.cumsum()
    u = np.random.rand(5000)
    choices = (u[..., np.newaxis] < c).argmax(axis=1)
    TEST_STR = "".join(np.array(["A", "T", "G", "C"])[choices])
    return TEST_STR

def test_convertKmerToLocation():
    """ Produce a correct insert location.
    """
    fn = kmer.create_kmer_loc_fn(7)

    assert fn("AAAAAA" + "A") == 0
    assert fn("AAAAAA" + "C") == 1
    assert fn("AAAAAA" + "G") == 2
    assert fn("AAAAAA" + "T") == 3

def test_encode():
    """ Produce a correct insert location
    """
    assert kmer.encode("AAAAAA" + "A") == (0, 7)
    assert kmer.encode("AAAAAA" + "C") == (1, 7)
    assert kmer.encode("AAAAAA" + "G") == (2, 7)
    assert kmer.encode("AAAAAA" + "T") == (3, 7)

def test_decode():
    """ It should be an inverse of encode function
    """
    for i in range(10):
        TEST_STR = _generate_sequence(l=8)
        assert kmer.decode(*kmer.encode(TEST_STR)) == TEST_STR

def test_guessFront():
    """ Guess nucleotide in front.
    """
    # Testing all edge case.
    q = {"AAAAAAAA": 0,
         "ATTTTTTT": 0,
         "CAAAAAAA": 1,
         "CTTTTTTT": 1,
         "GAAAAAAA": 2,
         "GTTTTTTT": 2,
         "TAAAAAAA": 3,
         "TTTTTTTT": 3}

    for s, a in q.items():
        h, l = kmer.encode(s)
        guess = np.digitize([h], kmer.guess_front(l))[0]
        assert guess == a


def testTrimFront():
    """ Trimming from front of kmer.
    """
    for i in ["A", "T", "G", "C"]:
        for j in ["AAAAAA", "TTTTTT"]:
            TEST_STR = i + j
            h, s = kmer.encode(TEST_STR)
            fTrim = kmer.trimFront(h, s)
            assert kmer.decode(fTrim, s - 1) == j
    # TEST_STRS = ["AAGACCAGATAC", "CAGACCAGATAC", "GAGACCAGATAC", "TAGACCAGATAC"]
    # for TEST_STR in TEST_STRS:
    #     h, s = kmer.encode(TEST_STR)
    #     fTrim = kmer.trimFront(h, s)
    #     assert kmer.decode(fTrim, s - 1) == "AGACCAGATAC"

    # Check for the edge case of
    TEST_STR = "G" + ("T" * 28)
    h, s = kmer.encode(TEST_STR)
    fTrim = kmer.trimFront(h, s)
    assert kmer.decode(fTrim, s - 1) == "T" * 28

def testTrimBack():
    """ Trimming from back of kmer.
    """
    TEST_STRS = ["TAGACCAGATAA", "TAGACCAGATAC", "TAGACCAGATAG", "TAGACCAGATAT"]
    for TEST_STR in TEST_STRS:
        h, s = kmer.encode(TEST_STR)
        bTrim = kmer.trimBack(np.array([h]))
        bStr = kmer.decode(bTrim[0], s - 1)
        assert bStr == "TAGACCAGATA"


def test_emit():
    """ Test the emit function.
    """
    TEST_STR = "ATGAATAGAGACTAGCATCTAGCTACGAT"
    f = list(kmer.kmer_hash_emit(TEST_STR, 5))
    s = list(kmer.kmer_hash_emit_rolling(TEST_STR, 5))
    assert f == s

    TEST_STR = "GGTTGGTGTGTGTGTGTGTGTTGAAA"
    f = list(kmer.kmer_hash_emit(TEST_STR, 4))
    s = list(kmer.kmer_hash_emit_rolling(TEST_STR, 4))
    assert f == s

def test_check_cut():
    """
    """
    for i in range(5):
        TEST_STR = _generate_sequence(l=500)
        KMER_B = 10
        KMER_S = KMER_B - 1
        bighashes = np.array(list(kmer.kmer_hash_emit(TEST_STR, KMER_B)))
        smallHashes = np.array(list(kmer.kmer_hash_emit(TEST_STR, KMER_S)))
        fHashes = kmer.trimFront(bighashes, KMER_B)
        bHashes = kmer.trimBack(bighashes)
        assert np.all(np.in1d(fHashes, smallHashes))
        assert np.all(np.in1d(bHashes, smallHashes))

def test_compare_encode():
    """ Test if encode is equal to create_kmer_loc_fn
    """
    test_fn1 = lambda x: kmer.encode(x)[0]
    test_fn2 = kmer.create_kmer_loc_fn(9)

    input1 = "ATGACGAGT"
    assert test_fn1(input1)  == test_fn2(input1)

    input2 = "TAAAAAAAA"
    assert test_fn1(input2)  == test_fn2(input2)

    input3 = "CAGAGCTAG"
    assert test_fn1(input3)  == test_fn2(input3)
