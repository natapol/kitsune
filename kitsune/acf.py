from . import kitsunejf as jf
import math
from tqdm import tqdm

def cal_acf(fsas, kmers, **karg):
    """ Calculate Average number of common features (ACF)
        acf(l) = sum((common feature of genome i and j at length l)/number of genome - 1)
    Args:
        fsas a list of genome file
        kmers a list of kmer to calculate

    Returns: dict(kmer: acf)

    """
    n = len(fsas)
    if n >= 2:
        result = {}
        for kmer in tqdm(kmers):
            keys_array = []
            for fsa in fsas:
                keys_array.append(set(jf.Kmercount(fsa, kmer, **karg).keys()))
            ccf = 0
            for pri_idx, key in enumerate(keys_array):
                for sec_idx in range(pri_idx + 1, n):
                    ccf += len(keys_array[pri_idx] & keys_array[sec_idx])
            result[kmer] = ccf/(n-1)
        return result
    else:
        raise

if __name__ == "__main__":
    print(cal_acf(['./examples/ASM170810v1_genomic.fna', './examples/ASM170810v1_genomic.fna'], [10,11]))
