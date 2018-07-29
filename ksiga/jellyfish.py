#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" use jellyfish to manipuate kmer
"""

import sys, os
import platform
import itertools
import subprocess
import numpy as np
from pathlib import Path


projectpath = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]

if platform.system() == 'Darwin':
    jellyfishpath = projectpath + "/bin/jellyfish-macosx"
elif platform.system() == 'Linux':
    jellyfishpath = projectpath + "/bin/jellyfish-linux"
else:
    pass 

def help(cmd = ''):
    return subprocess.getoutput('{} {} --help'.format(jellyfishpath, cmd))

def count(input, indexname, kmer):
    input_file = Path(input)
    if input_file.exists():
        subprocess.run('{0} count -m {1} -s {2} -o {3}_{1}.jf {4}'.format(jellyfishpath, kmer, '100M', output, input), shell=True)
    else:
        raise FileNotFoundError("{} is not found.".format(input))

def query(seq, indexname):
    file_name = "{}_{}.jf".format(indexname, len(seq))
    index_file = Path(file_name)
    if index_file.exists():
        return subprocess.getoutput('{} query {} {}'.format(jellyfishpath, index_file, seq)).split('\n')[0].split(' ')[1]
    else:
        raise FileNotFoundError("The jellyfish index file; {}; is not found.".format(file_name))

def queries(seq, *indexnames):
    counts = list()
    for idxn in indexnames:
        counts.append(query(seq, idxn))
    return np.asarray(counts, dtype=np.float64)

if __name__ == "__main__":
    # print(help())
    # print(count('{}/examples/S288C_reference_sequence_R64-2-1_20150113.fsa'.format(projectpath), 'S288C', 24))
    # print(query('S288C', 'AAAAAAAAAAAAAAAAAAAAAAAA'))
    idxnames = ['S288C'] * 28
    result = np.zeros(len(idxnames))
    for i in itertools.product('ATCG', repeat=24):
        result += queries(''.join(i), *idxnames)

    print(list(result))