#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Calculate and manipuate kmer
"""

import re
import itertools
import math
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
        seq (string): nucleotide sequence without N.
        ksize (int): size of kmer to calculate
        keyfn (func): hash function.

    Returns: TODO

    """
    for i in range(0, len(seq)-ksize+1):
        kmer = seq[i: i+ksize]
        #if not bool(NUCLEOTIDE_REGEX.search(kmer)):  # If it has no invalid character.
        yield keyfn(kmer)


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
    fn = lambda s: encode(s)[0]
    for seq in handle:
        rawStr = str(seq.seq).upper()
        rawSeqs = notAllowReg.split(rawStr)
        for rawSeq in rawSeqs:
            if len(rawSeq) < ksize:
                continue
            for hsh in kmer_hash_emit(rawSeq, ksize, fn):
                emits[locpointer] = hsh
                locpointer += 1
    
    indices, data = np.unique(emits[0:locpointer], return_counts=True)

    return  indices, data


# def kmer_count(seq, ksize, keyfn=id_fn):
#     """ Calculate kmer and store to a dictionary.

#     Args:
#         seq (string): nucleotide sequence
#         ksize (int): size of kmer to calculate
#         keyfn (func): key transformation

#     Returns: TODO

#     """
#     DICT = Counter()
#     for i in range(0, len(seq)-ksize+1):
#         kmer = seq[i: i+ksize]
#         if not bool(NUCLEOTIDE_REGEX.search(kmer)):  # If it has no invalid character.
#             DICT[keyfn(kmer)] += 1

#     return DICT

# #
# # Conveniant function
# #
# def kmer_count_fasta(f, ksize=13, keyfn=id_fn):
#     mainCounter = Counter()
#     handle = SeqIO.parse(f, "fasta")
#     for seq in handle:
#         seq = str(seq.seq)
#         mainCounter += kmer_count(seq, ksize, keyfn)
#     return mainCounter


# def build_csr_matrix_from_fasta(fh, ksize):
#     """ Build a csr row from fasta file.
#         Could be later build into a matrix.

#     Args:
#         f (str): fasta file to import

#     Returns: TODO

#     """

#     mainCounter = Counter()
#     handle = SeqIO.parse(fh, "fasta")
#     keyfn = create_kmer_loc_fn(ksize)

#     for seq in handle:
#         seq = str(seq.seq)
#         mainCounter += kmer_count(seq, ksize, keyfn)

#     # Initialize a sparse matrix.
#     rowNum = 1
#     colNum = kmer_location("A" * (ksize+1)) - kmer_location("A" * ksize)
#     shape = (rowNum, colNum)

#     tupCounter = ((k, v) for k, v in mainCounter.items())
#     indices, data = np.dstack(mainCounter.items())[0]  # Similar effect to zip(*item)
#     # Sort index so that we can efficiency search.
#     newIdx = np.argsort(indices)
#     indices = indices[newIdx]
#     data = data[newIdx]
#     return indices, data