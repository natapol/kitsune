"""
.. module:: kitsunejf
   :platform: Unix, MacOSX
   :synopsis: interface function to jellyfish

.. moduleauthor:: Natapol Pornputtapong <natapol.p@chula.ac.th>


"""
import os
import sys
import copy
import math
import shutil
import platform
import tempfile
import warnings
import subprocess
import collections

import numpy as np
from tqdm import tqdm
from scipy.stats import chi2
from scipy.spatial import distance
from scipy.special import loggamma


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

class NoInput(Exception):
    """Do X and return a list."""
    pass

if platform.system() == 'Darwin':
    jellyfishpath = os.path.join(__location__, 'modules', "jellyfish-macosx")
elif platform.system() == 'Linux':
    jellyfishpath = os.path.join(__location__, 'modules', "jellyfish-linux")
else:
    raise Exception("Windows are not supported.")

def mash(u, v, kmer):
    """Do X and return a list."""
    j = 1 - distance.jaccard(u, v)
    j = 1e-100 if j <= 0 else j
    return (-1/kmer) * math.log(2*j/(1+j))

def jsmash(u, v, kmer):
    """Do X and return a list."""
    j = 1 - distance.jensenshannon(u, v)
    j = 1e-100 if j <= 0 else j
    return (-1/kmer) * math.log(2*j/(1+j))


def nCr(n, r):
    """
    Calculate combinatorial using gamma funcion for huge number
    """
    logncr = loggamma(n+1) - loggamma(r+1) - loggamma(n-r+1)
    return math.exp(logncr)

def jaccarddistp(u, v):
    """
    Do X and return a list.
    """
    m = len(u)
    pu = u.sum()/m
    pv = v.sum()/m

    degenerate = False

    expectation = ((pu*pv)/(pu+pv-pu*pv))

    j_obs = np.logical_and(u, v).sum() / np.logical_or(u, v).sum() - expectation

    if(pu == 1 or pv == 1 or u.sum() == len(u) or v.sum() == len(v)):
        warnings.warn("One or both input vectors contain only 1's.", Warning)
        degenerate = True

    if(pu == 0 or pv == 0 or u.sum() == 0 or v.sum() == 0):
        warnings.warn("One or both input vectors contain only 0's", Warning)
        degenerate = True

    if degenerate:
        return 1.0

    # tan = jaccard_mca_rcpp(px,py,m,j.obs,accuracy)
    # pvalue = tan$pvalue

    #   pvalue <- switch(error.type, lower = pvalue, average = pvalue/tan$accuracy,
    #                upper = pvalue + 1 - tan$accuracy)

    # return(
    #     list(
    #     statistics = j.obs,
    #     pvalue = pvalue,
    #     expectation = expectation,
    #     accuracy = 1 - tan$accuracy,
    #     error.type = error.type
    #     )
    # )
    # Compute p-value using an asymptotic approximation

    q = [pu*pv, pu+pv-2*pu*pv]
    qq = q[0] + q[1]
    sigma = q[0] * q[1] * (1-q[0]) / (qq ** 3)
    norm = math.sqrt(m) * j_obs / math.sqrt(sigma)
    return chi2.pdf(norm*norm, 1)



DISTANCE_FUNCTION = {
    'braycurtis': distance.braycurtis,
    'canberra': distance.canberra,
    'chebyshev': distance.chebyshev,
    'cityblock': distance.cityblock,
    'correlation': distance.correlation,
    'cosine': distance.cosine,
    'dice': distance.dice,
    'euclidean': distance.euclidean,
    'hamming': distance.hamming,
    'jaccard': distance.jaccard,
    'kulsinski': distance.kulsinski,
    'rogerstanimoto': distance.rogerstanimoto,
    'russellrao': distance.russellrao,
    'sokalmichener': distance.sokalmichener,
    'sokalsneath': distance.sokalsneath,
    'sqeuclidean': distance.sqeuclidean,
    'yule': distance.yule,
    'jensenshannon': distance.jensenshannon,
    'mash': mash,
    'jsmash' : jsmash,
    'jaccarddistp': jaccarddistp
}

NUMERIC_DISTANCE = [
    distance.braycurtis,
    distance.canberra,
    distance.chebyshev,
    distance.cityblock,
    distance.correlation,
    distance.cosine,
    distance.euclidean,
    distance.sqeuclidean
]

BOOLEAN_DISTANCE = [
    distance.dice,
    distance.hamming,
    distance.jaccard,
    jaccarddistp,
    distance.kulsinski,
    distance.rogerstanimoto,
    distance.russellrao,
    distance.sokalmichener,
    distance.sokalsneath,
    distance.yule
]

PROB_DISTANCE =[
    distance.jensenshannon
]

class Kmercount(collections.Counter):
    """Do X and return a list."""
    def __init__(self, fsa, k_mer, **karg):

        self.kmer = k_mer

        if 'thread' not in karg:
            karg['thread'] = 1
        if 'lower' not in karg:
            karg['lower'] = 1
        if 'bchashsize' not in karg: #hashsize for jellyfish bc step
            karg['bchashsize'] = '1G'
        if 'hashsize' not in karg: #hashsize for jellyfish count step
            karg['hashsize'] = '100M'
        if 'canonical' not in karg or not karg['canonical']:
            canonical = ''
        else:
            canonical = '-C'

        if not os.path.isfile(fsa):
            raise NoInput('input is missing')

        filebasename = os.path.basename(fsa)

        with tempfile.TemporaryDirectory() as tmpdirname:

            if 'fast' in karg and karg['fast']:
                # for genome with one step
                dumpdata = subprocess.getoutput("""
                    {0} count {8} -m {1} -s {3} -t {4} {6} -o {5}.jf
                    {0} dump -c -L {7} {5}.jf
                    """.format(
                        jellyfishpath,
                        self.kmer,
                        karg['bchashsize'],
                        karg['hashsize'],
                        karg['thread'],
                        os.path.join(tmpdirname, filebasename),
                        fsa,
                        karg['lower'],
                        canonical
                    )
                )
            else:
                dumpdata = subprocess.getoutput("""
                    {0} bc {8} -m {1} -s {2} -t {4} -o {5}.bc {6}
                    {0} count {8} -m {1} -s {3} -t {4} --bc {5}.bc {6} -o {5}.jf
                    {0} dump -c -L {7} {5}.jf
                    """.format(
                        jellyfishpath,
                        self.kmer,
                        karg['bchashsize'],
                        karg['hashsize'],
                        karg['thread'],
                        os.path.join(tmpdirname, filebasename),
                        fsa,
                        karg['lower'],
                        canonical
                    )
                )

        datadict = {}
        for line in dumpdata.split('\n'):
            dat = line.rstrip().split(' ')
            datadict[dat[0]] = int(dat[1])

        super(Kmercount, self).__init__(datadict)

        # assign instance variable
        self.sum = sum(self.values())
        self.name = '.'.join(filebasename.split('.')[0:-1]).replace(' ', '_')

    def __repr__(self):
        return self.name

    def dist(self, other, dist_func):
        """Do X and return a list."""
        a, b = self.norm(other)
        if dist_func is mash:
            return dist_func(a.astype(bool), b.astype(bool), self.kmer)
        elif dist_func is jsmash:
            return dist_func(a.astype(float)/a.sum(), b.astype(float)/b.sum(), self.kmer)
        elif dist_func in BOOLEAN_DISTANCE:
            return dist_func(a.astype(bool), b.astype(bool))
        elif dist_func in PROB_DISTANCE:
            return dist_func(a.astype(float)/a.sum(), b.astype(float)/b.sum())
        else:
            return dist_func(a, b)

    def norm(self, other):
        """Do X and return a list."""
        mers = list(self.keys())
        mers.extend(list(other.keys()))
        mers = list(set(mers))
        a = []
        b = []
        for mer in mers:
            a.append(self[mer])
            b.append(other[mer])
        return np.array(a), np.array(b)

if __name__ == "__main__":
    pass
    # x = [0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1,
    #     1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0,
    #     0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]
    # y = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1,
    #     1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
    #     1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1]
    # print(jaccarddistp(np.array(x).astype(bool), np.array(y).astype(bool)))
    # print(prob_kmer_in_genome(16, 3e9))
    # for i in range(5, 15):
    #     genomea = Kmercount('../examples/S288C_reference_sequence_R64-2-1_20150113.fsa', i)
    #     genomeb = Kmercount('../examples/ASM170810v1_genomic.fna', i)
    #     print(genomea.dist(genomeb, mash))
    #     print(i, genomea.dist(genomeb, jaccarddistp))
    # for distfunc in NUMERIC_DISTANCE:
    #     print(genomea.dist(genomeb, distfunc))
    # for distfunc in BOOLEAN_DISTANCE:
    #     print(genomea.dist(genomeb, distfunc))
    # for distfunc in PROB_DISTANCE:
    #     print(genomea.dist(genomeb, distfunc))
