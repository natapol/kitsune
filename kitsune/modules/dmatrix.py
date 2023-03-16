import argparse as ap
import collections
import json
import sys
from operator import attrgetter

from tqdm import tqdm

from kitsune.modules import kitsunejf as jf


def read_params(args):
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog="kitsune (dmatrix)",
        description=(
            "Calculate a distance matrix"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--filenames", nargs="*", help="Genome files in fasta format")
    p.add_argument("--fast", action="store_true", help="Jellyfish one-pass calculation (faster)")
    p.add_argument("--canonical", action="store_true", help="Jellyfish count only canonical mer")
    p.add_argument("-k", "--kmer", required=True, type=int)
    p.add_argument("-i", "--input", type=str, help="List of genome files in txt")
    p.add_argument("-o", "--output", type=str, help="Output filename")
    p.add_argument("-t", "--thread", type=int, default=1)
    p.add_argument("--transformed", action="store_true")
    p.add_argument(
        "-d", 
        "--distance", 
        default="cosine", 
        choices=["braycurtis", "canberra", "jsmash", "chebyshev", "cityblock", 
                 "correlation", "cosine", "dice", "euclidean", "hamming", 
                 "jaccard", "kulsinsk", "matching", "rogerstanimoto", "russellrao", 
                 "sokalmichener", "sokalsneath", "sqeuclidean", "yule", "mash", "jaccarddistp"])
    p.add_argument("-f", "--format", default="phylip")
    return p.parse_args(args)


def cal_matrix(fsas, k_mer, dist_func, **karg):
    result = collections.OrderedDict()
    counters = list()
    n = len(fsas)
    
    # create kmercounter
    for fsa in tqdm(fsas, desc="k-mer counting"):
        counters.append(jf.Kmercount(fsa, k_mer, **karg))

    sorted(counters, key=attrgetter('name'))

    for counter in counters:
        result[counter.name] = [0] * n

    for pri_idx, pri_counter in tqdm(enumerate(counters), desc="distant calculation"):
        for sec_idx in range(pri_idx + 1, n):
            sec_counter = counters[sec_idx]
            dist = pri_counter.dist(sec_counter, dist_func, transform=karg['transformed'])
            result[pri_counter.name][sec_idx] = dist
            result[sec_counter.name][pri_idx] = dist
    
    return result


def run(args):
    """
    Calculate a distance matrix
    """

    # Load command line parameters
    args = read_params(args)

    args.distance_func = jf.DISTANCE_FUNCTION[args.distance]

    if args.input:
        genomeList = open(args.input, "r").read().strip().split("\n")
        outdata = cal_matrix(genomeList, args.kmer, jf.DISTANCE_FUNCTION[args.distance], **vars(args))
    
    else:
        outdata = cal_matrix(args.filenames, args.kmer, jf.DISTANCE_FUNCTION[args.distance], **vars(args))
    
    if args.format == "phylip":
        outdata = str(len(outdata)) + '\n' + '\n'.join([k + '\t' + '\t'.join([str(x) for x in v]) for k, v in outdata.items()])
    
    else:
        outdata = json.dumps(outdata, indent=2)

    print(outdata, file=open(args.output, "w+") if args.output else None)


if __name__ == "__main__":
    run(sys.argv)
