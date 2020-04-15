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
kopt <filenames>                  Compute optimal kmer and generate histograms .
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
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-ke", "--kend", required=True, type=int, help="last k-mer")
    parser.add_argument("-kf", "--kfrom", default=4, type=int, help="Calculate from k-mer")
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
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
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
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmers", nargs="+", required=True, type=int)
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
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-k", "--kmer", required=True, type=int)
    parser.add_argument("-i", "--input", type=str, help="list of genome files in txt")
    parser.add_argument("-o", "--output", type=str, help="output filename")
    parser.add_argument("-t", "--thread", type=int, default=1)
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

##################################### New Additions ###################################################

def filter(n,m): #### A filter to select kmers that meet a criteria ########################
    c = Counter(n)
    cr = set(n)
    result = []
    for i in cr:
        if c[i] >= (m * 0.5):
            result.append(i)
    return result

def kopt_cre(l,genomeList,args): ###################### Kmer optimal cre #######################
    from . import kopt
    genome = genomeList[l]
    outdata = kopt.cal_cre(genome, **vars(args))
    outdata = sorted(outdata.items(), key=itemgetter(0))
    return outdata

def kopt_acf(p,args):           ##################### kmer optimal acf ########################
    from . import kopt
    outdata_acf = kopt.cal_acf(p, **vars(args))
    outdata_acf = sorted(outdata_acf.items(), key=itemgetter(0))
    return outdata_acf


def kopt_ofc(genomeList,kmer, args):
    from . import kopt
    outdata_ofc = kopt.cal_ofc(genomeList,kmer, **vars(args))
    outdata_ofc = sorted(outdata_ofc.items(), key=itemgetter(0))
    return outdata_ofc

def chunk(it, size):
    it = iter(it)
    return iter(lambda: tuple(islice(it, size)), ())

def final_results(acf_res, outdata_cre, outdata_ofc): #################### some calculations for acf & cre #########
    outdata_ac, outdata_cr, res, result, outdata_of = ([] for i in range(5))
    outdata_of= '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata_ofc])
    ofc_max = max([x[1] for x in outdata_ofc])
    kmax = [x for x,y in outdata_ofc if y==ofc_max]
    ofc_possible_kmers = [x[0] for x in outdata_ofc if x[0] >= kmax[0]]
    outdata_ac= '\n'.join(['\t'.join([str(x) for x in data]) for data in acf_res])
    k_acf = [x[0] for x in acf_res]
    acf = [x[1] for x in acf_res]
    acf_max = max(acf)
    maxacf = 0.10*float(acf_max)
    acf_kmax = [x[0] for x in acf_res if x[1] <= maxacf and x[1] >= 1.0]
    outdata_cr= '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata_cre])
    cre = [x[1] for x in outdata_cre]
    if len(cre) != 0:
        maxval = 0.1*float(max(cre))
        result = [x[0] for x in outdata_cre if x[1] <= maxval and x[1] >= 1e-4]
    else:
        pass
    return acf_kmax, result, outdata_cr, outdata_ac, outdata_of, ofc_possible_kmers;


def kmer_optimal(args):               ################## kmer_optimal #################################
    from . import kopt
    """ Find the recommended optimal kmer using acf,cre and ofc for a given set of genomes

    Args:
        args (TODO): TODO

    Returns: TODO

    """
    desc = "Find the recommended choice of kmer[s] for a given set of genomes within a given interval"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filenames", nargs="*",type=str, help="genome files in fasta format") #nargs"*" is added 
    parser.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    parser.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer (use for raw read count)")
    parser.add_argument("-ks", "--ksmall", default=5, type=int, help="smallest k-mer length to consider")
    parser.add_argument("-kl", "--klarge", required=True, type=int, help="largest k-mer length to consider")
    parser.add_argument("-i", "--input", type=str, help="path to list of genome files in a text format") # input introduced
    parser.add_argument("-o", "--output", type=str, help="output filename")
    args = parser.parse_args(args)

    if args.input:
        genomeList = open(args.input, "r").read().strip().split("\n")
        combs = combinations(genomeList,2)
        genLen = len(genomeList)
        combLen = float(len(list(combinations(genomeList,2))))
    
        ######################## multiprocessing for cre ####################################
        outdata_cre = []
        acf_res = []
        pool = mp.Pool()
        ##############changes        

        acf_res = [pool.apply_async(kopt_acf, args=(p,args)) for p in combs]
        outdata_cre = [pool.apply_async(kopt_cre, args=(l,genomeList,args)) for l in range(genLen)]
           
        ######i################# observed feature occurrence #################################

        for kmer in range(args.ksmall, args.klarge+1):
            outdata_ofc.append(pool.apply_async(kopt_ofc, args=(genomeList, kmer,args)))
        pool.close()
        pool.join()
    outdata_ofc = list(chain(*[result.get() for result in outdata_ofc]))
    outdata_cre = list(chain(*[result.get() for result in outdata_cre]))
    acf_res = list(chain(*[result.get() for result in acf_res]))
    outdata_of, outdata_cr, outdata_ac = ([] for i in range(3)) 
    outdata_of= '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata_ofc])
    outdata_ac= '\n'.join(['\t'.join([str(x) for x in data]) for data in acf_res])
    outdata_cr= '\n'.join(['\t'.join([str(x) for x in data]) for data in outdata_cre])
    final = final_results(acf_res, outdata_cre, outdata_ofc)
    result2 = filter(final[1], genLen)
    acf_kmax2 = filter(final[0], combLen)
    opt = list(set(result2) & set(acf_kmax2) & set(final[5]))
     
    if args.output is not None and len(opt)>0:
        with open(args.output, 'w') as ofhandle:
            print("The recommended choice of kmer: ", min(opt), file=ofhandle)
            print("OFC possible kmers= :",final[5], file=ofhandle)
            print("CRE possible kmers= :",list(set(result2)), file=ofhandle)
            print("ACF possible kmers= :",list(set(acf_kmax2)), file = ofhandle)
    elif args.output is not None and len(opt)==0:
        with open(args.output, 'w') as ofhandle:
            print("The recommended choice of kmer does not lie within the given interval", file=ofhandle)
            print("OFC possible kmers :",final[5], file=ofhandle)
            print("CRE possible kmers= :",list(set(result2)), file=ofhandle)
            print("ACF possible kmers= :",list(set(acf_kmax2)), file = ofhandle)
    elif len(opt) == 0:
        print("The recommended choice of kmer does not lie within the given interval")
        print("OFC possible kmers :",final[5])
        print("CRE possible kmers= :",list(set(result2)))
        print("ACF possible kmers= :",list(set(acf_kmax2)))
    else:
        print("The recommended choice of kmer: ",min(opt))
   

if __name__ == "__main__":
    main()
