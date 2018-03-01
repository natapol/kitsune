#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and manipuate kmer
"""

import re
import itertools
import math
from collections import deque
from collections import Counter

import numpy as np
import Bio.SeqIO as SeqIO

# Global variable
KMER_ARR = ["A", "C", "G", "T"]
KSIZE_LIMIT = 31

# Regex to clean sequence data
NUCLEOTIDE_REGEX = re.compile("[^ATGC]")

# Function to collapse kmer

id_fn = lambda x: x


def kmer_location(kmer):
    """ Calculate kmer location to store in array. Note that AAAA does not end up in 0 since with this encoding scheme
        leave some space for lower kmer.
        NOTE: Not compatible with decode.
    """
    encMap = {"A": 1, "C": 2, "G": 3, "T": 4}

    code = 0
    for ch in kmer:
        code *= 4
        code += encMap[ch]

    return code


def encode(k):
    """
    """
    encoding_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map[ch]
    return code, len(k)


def decode(code, length):
    """ Reverse
    """
    decoding_lst = ["A", "C", "G", "T"]
    ret = ''
    for _ in range(length):
        index = code & 3
        code >>= 2
        ret = decoding_lst[index] + ret
    return ret


def guess_front(ksize):
    """
    """
    # Guess a front of character, return a result that we could use np.digitize
    # bins = guess_front(8)
    # np.digitize([hash], bins) -> 0 = A, 1 = C, 2 = G, 3 = T

    A_lim = encode("A" + ("T" * (ksize - 1)))[0]
    C_lim = int((A_lim * 2) + 1)
    G_lim = int((A_lim * 3) + 2)
    return np.array([A_lim + 1, C_lim + 1, G_lim + 1])


def trimFront(khashs, ksize):
    """ Calculate a hash of k-mer if the last character is trimmed off.

    Args:
        khashs (np.array[Int])

    Returns:
        khashs (np.array[Int])
    """
    bins = guess_front(ksize)
    if ksize < 27:
        frontCharHash = np.digitize(khashs, bins)
    else:  # digitize have problem with a very large number.
        frontCharHash = np.searchsorted(bins, khashs, side='right')
    fHash = khashs - (frontCharHash * (4 ** (ksize - 1)))
    return fHash


def trimBack(khashs):
    """ Calculate a hash of k-mer if the last character is trimmed off.

    Args:
        khashs (np.array[Int])

    Returns:
        khashs (np.array[Int])
    """
    return khashs // 4


def generateMers(size=4):
    """Return a list of mers
    """
    kmerList = KMER_ARR
    it = itertools.product(kmerList, repeat=size)
    for mers in it:
        yield "".join(mers)


def create_kmer_loc_fn(size):
    """ Hash location of kmer for specific size.
    NOTE: This is pretty much similar to encode. May refactor later.
    """

    offset = kmer_location("A" * size)

    def wrapped(seq):
        return kmer_location(seq) - offset

    return wrapped


def kmer_hash_emit(seq, ksize, keyfn=encode):
    """ Calculate kmer and emit them.

    Args:
        seq (string): Cleaned, uppercase nucleotide sequence.
        ksize (int): size of kmer to calculate
        keyfn (func): hash function.

    Returns: Yield hash kmer in each.

    """
    for i in range(0, len(seq) - ksize + 1):
        kmer = seq[i: i + ksize]
        yield keyfn(kmer)[0]


def kmer_hash_emit_rolling(seq, ksize):
    """ Calculate kmer's hash and emit it.
    This one is specialize and optimize by using rolling hashing method.

    Args:
        seq (string): nucleotide sequence without N.
        ksize (int): size of kmer to calculate
        keyfn (func): hash function.

    Returns: Yield hash kmer in each.

    """
    encoding_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    # Initialize a queue
    queue = deque(seq[0:ksize], ksize)
    hash_val = encode(queue)[0]
    yield hash_val

    for char in seq[ksize:]:
        # Delete front and add back
        # update hash
        hash_val -= encoding_map[queue[0]] * (4 ** (ksize - 1))
        hash_val *= 4  # Shift
        queue.append(char)
        hash_val += encoding_map[queue[-1]]
        yield hash_val