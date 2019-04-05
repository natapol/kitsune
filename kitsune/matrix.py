import kitsunejf as jf
import math
import collections
import numpy as np
from tqdm import tqdm
from operator import attrgetter

def cal_matrix(fsas, kmer, dist_func=jf.DISTANCE_FUNCTION['cosine'], **karg):
    result = collections.OrderedDict()
    counters = []
    n = len(fsas)
    # create kmercounter
    for fsa in tqdm(fsas):
        counters.append(jf.Kmercount(fsa, **karg))

    sorted(counters, key=attrgetter('name'))

    for counter in counters:
        result[counter.name] = [0] * n

    for pri_idx, pri_counter in enumerate(counters):
        for sec_idx in range(pri_idx + 1, n):
            sec_counter = counters[sec_idx]
            dist = pri_counter.dist(sec_counter, dist_func)
            result[pri_counter.name][sec_idx] = dist
            result[sec_counter.name][pri_idx] = dist
    return result

if __name__ == "__main__":
    print(cal_matrix(['./examples/S288C_reference_sequence_R64-2-1_20150113.fsa', './examples/ASM170810v1_genomic.fna'], 25))