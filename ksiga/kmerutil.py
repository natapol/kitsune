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

# Regex to clean sequence data
NUCLEOTIDE_REGEX = re.compile("[^ATGC]")

# Function to collapse kmer
id_fn = lambda x: x

def kmer_location(kmer):
    """ Calculate kmer location to store in array. Note that AAAA does not end up in 0 since with this encoding scheme
        leave some space for lower kmer.
        NOTE: Not compatible with decode.
    """
    encMap = {"A":1, "C":2, "G":3, "T":4}

    code = 0
    for ch in kmer:
        code *= 4
        code += encMap[ch]

    return code


def encode(k):
    """
    """
    encoding_map = {"A":0, "C":1, "G":2, "T":3}
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map[ch]
    return code, len(k)


def decode(code, length):
    """ Reverse
    """
    decoding_lst = KMER_ARR
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

    A_lim = encode("A" + ("T" * (ksize-1)))[0]
    C_lim = int((A_lim * 2) + 1)
    G_lim = int((A_lim * 3) + 2)
    return np.array([A_lim + 1, C_lim + 1, G_lim + 1])

def trimFront(khashs, ksize):
    """ Trim k-mer from front

    Args:
        khashs (np.array)
    """
    bins = guess_front(ksize)
    frontCharHash = np.digitize(khashs, bins)
    fHash = khashs - (frontCharHash * (4 ** (ksize - 1)))
    return fHash

def trimBack(khashs):
    """ Trim k-mer from back

    Args:
        khashs (np.array)
    """
    m = khashs % 4  # Determine what is in the back.
    bHashs = ((khashs - m)/4).astype(np.int64)
    return bHashs

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
    for i in range(0, len(seq)-ksize+1):
        kmer = seq[i: i+ksize]
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
    encoding_map = {"A":0, "C":1, "G":2, "T":3}
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


def build_csr_matrix_from_fasta(fh, ksize):
    # Check total number of base
    handleCount = SeqIO.parse(fh, "fasta")
    notAllowReg = re.compile("[^ATCG]")
    count = 1
    for seq in handleCount:
        count += len(seq)
    handleCount.close()
    emits = np.zeros(count) # At least should be THIS much.
    
    # Initialize iterator again, now to count
    fh.seek(0,0)
    handle = SeqIO.parse(fh, "fasta")
    locpointer = 0
    for seq in handle:
        rawStr = str(seq.seq).upper()
        rawSeqs = notAllowReg.split(rawStr)
        for rawSeq in rawSeqs:
            if len(rawSeq) < ksize:
                continue
            for hsh in kmer_hash_emit_rolling(rawSeq, ksize):
                emits[locpointer] = hsh
                locpointer += 1
    
    indices, data = np.unique(emits[0:locpointer], return_counts=True)

    return  indices, data