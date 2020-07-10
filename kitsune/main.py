#!/usr/bin/env python
"""
.. module:: main
   :platform: Unix, MacOSX
op involving two conditions python
   :synopsis: Main entry of the script.

.. moduleauthor:: Natapol Pornputtapong <natapol.p@chula.ac.th>


"""
import argparse
import sys
import os
import gzip
import datetime
import json
import collections
import numpy as np
from operator import itemgetter, attrgetter
from itertools import combinations, chain, islice
from collections import Counter
import multiprocessing as mp
import time

# import internal library
from . import kitsunejf as jf
from tqdm import *


USAGE = """
"""

def main():
    commands = {
                "cre": cumulative_relative_entropy,
                "acf": average_common_feature,
                "ofc": observe_feature_occurrence,
                "kopt": kmer_optimal,
                "dmatrix": generate_distance_matrix
               }

    parser = argparse.ArgumentParser(description="Compute k-mer signature of genome",
                                     usage="""kitsune <command> [<args>]

Commands can be:
cre <filename>                    Compute cumulative relative entropy.
acf <filenames>                   Compute average number of common feature between signatures.
ofc <filenames>                   Compute observed feature frequencies.
kopt <filenames>                  Compute recommended choice (optimal) of kmer within a given kmer interval for a set of genomes using the cre, acf and ofc.
dmatrix <filenames>               Compute distance matrix.
""")
    parser.add_argument('command')
    args = parser.parse_args(sys.argv[1:2])
    if args.command not in commands:
        parser.print_help()
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])

def cumulative_relative_entropy(args):
    from . import cre
    """ Calculate optimal k-mer through CRE value.

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate k-mer from cumulative relative entropy of all genomes"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", type=str, help="a genome file in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    parser.add_argument("-ke", "--kend", required=True, type=int, help="last k-mer")
    parser.add_argument("-kf", "--kfrom", default=4, type=int, help="Calculate from k-mer")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    outdata = cre.cal_cre(args.filename, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)


def average_common_feature(args):
    from . import acf
    """ Calculate an average number of common feature pairwise
        between one genome against others

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate average number of common feature"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", type=str, help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int, help="have to state before")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    outdata = acf.cal_acf(args.filenames, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

def observe_feature_occurrence(args):
    from . import ofc
    """ Calculate an observe feature frequency

    Args:
        args (TODO): TODO

    Returns: TODO

    """

    desc = "Calculate observe feature occurrence"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="+", help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    outdata = ofc.cal_ofc(args.filenames, args.kmers, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    outdata = '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata])
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

def generate_distance_matrix(args):
    from . import matrix
    """Generate distance matrix base on k-mer

    The output will

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Calculate a distance matrix"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="*", help="genome files in fasta format")
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    parser.add_argument("-k", "--kmer", required=True, type=int)
    parser.add_argument("-i", "--input", type=str, help="list of genome files in txt")
    parser.add_argument("-o", "--output", type=str, help="output filename")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("--transformed", action="store_true")
    parser.add_argument("-d", "--distance", default="cosine", help="braycurtis, canberra, jsmash, chebyshev, cityblock, correlation, cosine (default), dice, euclidean, hamming, jaccard, kulsinsk, matching, rogerstanimoto, russellrao, sokalmichener, sokalsneath, sqeuclidean, yule, mash, jaccarddistp")
    parser.add_argument("-f", "--format", default="phylip")
    args = parser.parse_args(args)
    args.distance_func = jf.DISTANCE_FUNCTION[args.distance]
    if args.input:
        genomeList = open(args.input, "r").read().strip().split("\n")
        outdata = matrix.cal_matrix(genomeList, args.kmer, jf.DISTANCE_FUNCTION[args.distance], **vars(args))
    else:
        outdata = matrix.cal_matrix(args.filenames, args.kmer, jf.DISTANCE_FUNCTION[args.distance], **vars(args))
    if args.format == "phylip":
        outdata = str(len(outdata)) + '\n' + '\n'.join([k + '\t' + '\t'.join([str(x) for x in v]) for k, v in outdata.items()])
    else:
        outdata = json.dumps(outdata, indent=2)
        
    if args.output is not None:
        with open(args.output, 'w') as ofhandle:
            ofhandle.write(outdata)
    else:
        print(outdata)

def filter(n,m): 
    c = Counter(n)
    cr = set(n)
    result = []
    for i in cr:
        if c[i] >= int(m * 0.5):
            result.append(i)
    return result

def kopt_cre(l, kl, genomeList, args):
    from . import kopt
    genome = genomeList[l]
    outdata = kopt.cal_cre(genome, kl, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    return outdata

def kopt_acf(p,args): 
    from . import kopt
    outdata_acf = kopt.cal_acf(p, **vars(args))
    outdata_acf = sorted(outdata_acf.items(), key=itemgetter(0))
    return outdata_acf

def kopt_ofc(genomeList,kmer, args):
    from . import kopt
    outdata_ofc = kopt.cal_ofc(genomeList,kmer, **vars(args))
    outdata_ofc = sorted(outdata_ofc.items(), key=itemgetter(0))
    return outdata_ofc

def Range(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def valid_file(x):
    base, ext = os.path.splitext(x)
    if ext.lower() not in ('.csv', '.tab','.txt',''):
        raise argparse.ArgumentTypeError('File must have a csv, txt, tab extension or no extension')
    return x

def kmer_optimal(args):
    from . import kopt
    """ Find the recommended optimal kmer using acf,cre and ofc for a given set of genomes

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Example: kitsune kopt genomeList.txt -kl 15 --canonical --fast -t 4 -o out.txt"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", type=valid_file, help="A file that list the path to all genomes(fasta format) with extension as (.txt,.csv,.tab) or no extension")
    parser.add_argument("--fast", action='store_true', help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action='store_true', help="Jellyfish count only canonical mer")
    parser.add_argument("-kl", "--klarge", required=True, type=int, help="largest k-mer length to consider, note: the smallest kmer length is 4")
    parser.add_argument("-o", "--output", type=str, help="output filename")
    parser.add_argument("--closely_related", action='store_true', help="For closely related set of genomes, use this option" )
    parser.add_argument("-x", "--cre_cutoff", type=Range, default= 0.1, help="cutoff to use in selecting kmers whose cre's are <= (cutoff * max(cre)), Default = 10 percent, ie x=0.1")
    parser.add_argument("-y", "--acf_cutoff", type=Range, default= 0.1, help="cutoff to use in selecting kmers whose acf's are >= (cutoff * max(acf)), Default = 10 percent, ie y=0.1")
    parser.add_argument("-t", "--thread", type=int, default=1, help="Number of threads (integer)")
    args = parser.parse_args(args)

    genomeList = open(args.filenames, "r").read().strip().split("\n")
    genLen = len(genomeList)

    if args.filenames and genLen >=2:
        combs = combinations(genomeList,2)
        combLen = float(len(list(combinations(genomeList,2))))
        xcutoff = args.cre_cutoff
        ycutoff = args.acf_cutoff
        outdata_cre = []
        acf_res = []
        pool = mp.Pool()

        if args.closely_related:
            kl = args.klarge
            acf_kmax2 = range(4, kl+1)
            pool1 = mp.Pool()
            for l in range(genLen):
                outdata_cre.append(pool1.apply_async(kopt_cre, args=(l, kl, genomeList,args)))
            outdata_cre = list(chain(*[result.get() for result in outdata_cre]))
            cre = [x[1] for x in outdata_cre]
            if len(cre) != 0:
                maxval = xcutoff * float(max(cre))
            else:
                maxval = 0
            result = [x[0] for x in outdata_cre if x[1] < maxval]
            result2 = filter(result, genLen * 0)

            pool1.close()
            pool1.join()
            time.sleep(10)

            pool2 = mp.Pool()
            outdata_ofc = []
            if len(result2) != 0:
                ks = min(result2)
            else:
                ks = 4
            for kmer in range(ks, kl+1):
                outdata_ofc.append(pool2.apply_async(kopt_ofc, args=(genomeList, kmer,args)))
            outdata_ofc = list(chain(*[result.get() for result in outdata_ofc]))
            ofc = [x[1] for x in outdata_ofc]
            if len(ofc) != 0:
                ofc_max = max(ofc)
            else:
                ofc_max = 1e7

            kmax = [x for x,y in outdata_ofc if y==ofc_max]
            ofc_possible_kmers = [x[0] for x in outdata_ofc if x[0] >= kmax[0]]
            pool2.close()
            pool2.join()

        else:
            for p in combs:
                acf_res.append(pool.apply_async(kopt_acf, args=(p,args)))
            acf_res = list(chain(*[result.get() for result in acf_res]))
            k_acf = [x[0] for x in acf_res]
            acf = [x[1] for x in acf_res]
            maxacf = ycutoff * float(max(acf))
            acf_kmax = [x[0] for x in acf_res if x[1] >= maxacf]
            acf_kmax2 = filter(acf_kmax, combLen * 0 )

            if len(acf_kmax2) == 0:
                kl = args.klarge
            else:
                kl =max(acf_kmax2)

            pool.close()
            pool.join()
            time.sleep(10)

            pool1 = mp.Pool()
            for l in range(genLen):
                outdata_cre.append(pool1.apply_async(kopt_cre, args=(l, kl, genomeList,args)))
            outdata_cre = list(chain(*[result.get() for result in outdata_cre]))
            cre = [x[1] for x in outdata_cre]
            if len(cre) != 0:
                maxval = xcutoff * float(max(cre))
            else:
                maxval = 0
            result = [x[0] for x in outdata_cre if x[1] < maxval]
            result2 = filter(result, genLen * 0)

            pool1.close()
            pool1.join()
            time.sleep(10)
            
            pool2 = mp.Pool()
            outdata_ofc = []
            if len(result2) != 0:
                ks = min(result2)
            else:
                ks = 4
            for kmer in range(ks, kl+1):
                outdata_ofc.append(pool2.apply_async(kopt_ofc, args=(genomeList, kmer,args)))
            outdata_ofc = list(chain(*[result.get() for result in outdata_ofc]))
            ofc = [x[1] for x in outdata_ofc]
            if len(ofc) != 0:
                ofc_max = max(ofc)
            else:
                ofc_max = 1e7

            kmax = [x for x,y in outdata_ofc if y==ofc_max]
            ofc_possible_kmers = [x[0] for x in outdata_ofc if x[0] >= kmax[0]]

            pool2.close()
            pool2.join()

    else:
        return print('list of filenames must be more than 1')

    opt = list(set(result2) & set(acf_kmax2) & set(ofc_possible_kmers))

    if args.output is not None and len(opt)>0:
        with open(args.output, 'w') as ofhandle:
            print("The recommended choice of kmer is: ", min(opt), file=ofhandle)
    elif args.output is not None and len(opt)==0:
        with open(args.output, 'w') as ofhandle:
            print("The recommended choice of kmer does not lie within the given interval from 4 to the largest kmer length you selected", file=ofhandle)
    elif len(opt) == 0:
        print("The recommended choice of kmer does not lie within the given interval from 4 to the largest kmer length you selected")

    else:
        print("The recommended choice of kmer is: ",min(opt))
            
if __name__ == "__main__":
    main()
