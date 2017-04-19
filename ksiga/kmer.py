#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and manipuate kmer
"""

import re
import itertools
from collections import Counter

import numpy as np
import Bio.SeqIO as SeqIO

# Global variable
KMER_ARR = ["A", "C", "G", "T"]
KMER_DIC = {"A":1, "C":2, "G":3, "T":4}

# Conveniant function
REGEX = re.compile("[^ATGC]")
id_fn = lambda x: x


def kmer_location(kmer):
    """ Calculate kmer location to store in array
    """
    encMap = KMER_DIC

    code = 0
    for ch in kmer:
        code *= 4
        code += encMap[ch]

    return code


def encode(k):
    encoding_map = KMER_DIC
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map[ch]
    return code, len(k)


def decode(code, length):
    decoding_lst = KMER_ARR
    ret = ''
    for _ in range(length):
        index = code & 3
        code >>= 2
        ret = decoding_lst[index] + ret
    return ret


def generateMers(size=4):
    """Return a list of mers
    """
    kmerList = KMER_ARR
    it = itertools.product(kmerList, repeat=size)
    for mers in it:
        yield "".join(mers)


def create_kmer_loc_fn(size):
    """ Hash location of kmer for specific size.
    """

    offset = kmer_location("A" * size)
    def wrapped(seq):
        return kmer_location(seq) - offset

    return wrapped

        
def kmer_count(seq, ksize, keyfn=id_fn):
    """ Calculate kmer and location to store

    Args:
        arg1 (TODO): TODO

    Returns: TODO

    """
    DICT = Counter()
    for i in range(0, len(seq)-ksize+1):
        kmer = seq[i: i+ksize]
        if not bool(REGEX.search(kmer)):
            DICT[keyfn(kmer)] += 1

    return DICT


#
# Conveniant function
#


def kmer_count_fasta(f, ksize=13, keyfn=id_fn):
    mainCounter = Counter()
    handle = SeqIO.parse(f, "fasta")
    for seq in handle:
        seq = str(seq.seq)
        mainCounter += kmer_count(seq, ksize, keyfn)
    return mainCounter


def build_csr_matrix_from_fasta(fh, ksize):
    """ Build a `row` of csr matrix from fasta file.
        Could be later build into a matrix.

    Args:
        f (str): fasta file to import

    Returns: TODO

    """

    mainCounter = Counter()
    handle = SeqIO.parse(fh, "fasta")
    keyfn = create_kmer_loc_fn(ksize)

    for seq in handle:
        seq = str(seq.seq)
        mainCounter += kmer_count(seq, ksize, keyfn)

    # From main counter to csr_matrix
    rowNum = 1
    colNum = kmer_location("A" * (ksize+1)) - kmer_location("A" * ksize)
    shape = (rowNum, colNum)

    tupCounter = ((k, v) for k, v in mainCounter.items())
    indices, data = np.dstack(mainCounter.items())[0]  # Similar effect to zip(*item)
    # We will exploit the fact that the indices is sorted in later search.
    newIdx = np.argsort(indices)
    indices = indices[newIdx]
    data = data[newIdx]
    return indices, data
