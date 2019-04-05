import os
import sys
import copy
import shutil
import platform
import tempfile
import subprocess
import collections
import numpy as np
from tqdm import tqdm
from scipy.spatial import distance

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

class NoInput(Exception):
    """Do X and return a list."""
    pass

if platform.system() == 'Darwin':
    jellyfishpath = os.path.join(__location__, 'modules', "jellyfish-macosx")
elif platform.system() == 'Linux':
    jellyfishpath = os.path.join(__location__, 'modules', "jellyfish-linux")
else:
    raise

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
    'matching': distance.matching,
    'rogerstanimoto': distance.rogerstanimoto,
    'russellrao': distance.russellrao,
    'sokalmichener': distance.sokalmichener,
    'sokalsneath': distance.sokalsneath,
    'sqeuclidean': distance.sqeuclidean,
    'yule': distance.yule
}

class Kmercount(collections.Counter):
    """Do X and return a list."""
    def __init__(self, fsa, kmer = 21, **karg):
        if 'thread' not in karg:
            karg['thread'] = 1
        if 'lower' not in karg:
            karg['lower'] = 1
        if 'bchashsize' not in karg: #hashsize for jellyfish bc step
            karg['bchashsize'] = '1G'
        if 'hashsize' not in karg: #hashsize for jellyfish count step
            karg['hashsize'] = '100M'
        if 'canonical' not in karg:
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
                        kmer,
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
                        kmer,
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

    def dist(self, other, dist_func = distance.cosine):
        """Do X and return a list."""
        a, b = self.norm(other)
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
    genomea = Kmercount('../examples/S288C_reference_sequence_R64-2-1_20150113.fsa')
    genomeb = Kmercount('../examples/ASM170810v1_genomic.fna')
    print(genomea.dist(genomeb))