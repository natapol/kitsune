#!/usr/bin/env python3
"""
Optimal kmer size selection for a set of genomes using Average number of Common Features (ACF),
Cumulative Relative Entropy (CRE), and Observed Common Features (OCF)
"""

import argparse as ap
import collections
import errno
import math
import multiprocessing as mp
import os
import sys
import time
from functools import partial
from itertools import chain, combinations
from operator import itemgetter
from typing import Dict, List, Optional, Tuple

import numpy as np  # type: ignore
import tqdm  # type: ignore

from kitsune.modules import kitsunejf as jf


def read_params(args):
    """
    Read and test input arguments

    :return:    The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog="kitsune (kopt)",
        description=(
            "Optimal kmer size selection for a set of genomes using Average number of Common Features (ACF), "
            "Cumulative Relative Entropy (CRE), and Observed Common Features (OCF). "
            "Example: kitsune kopt --filenames genomeList.txt --k-min 4 --k-max 12 --canonical --fast"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "--acf-cutoff",
        type=float,
        default=0.1,
        dest="acf_cutoff",
        help="Cutoff to use in selecting kmers whose ACFs are >= (cutoff * max(ACF))"
    )
    p.add_argument(
        "--canonical",
        action="store_true",
        default=False,
        help="Jellyfish count only canonical kmers"
    )
    p.add_argument(
        "--closely-related",
        action="store_true",
        default=False,
        dest="closely_related",
        help="Use in case of closely related genomes"
    )
    p.add_argument(
        "--cre-cutoff",
        type=float,
        default=0.1,
        dest="cre_cutoff",
        help="Cutoff to use in selecting kmers whose CREs are <= (cutoff * max(CRE))"
    )
    p.add_argument(
        "--fast",
        action="store_true",
        default=False,
        help="Jellyfish one-pass calculation (faster)"
    )
    p.add_argument(
        "--filenames",
        type=os.path.abspath,
        required=True,
        help=(
            "Path to the file with the list of genome files paths. "
            "There should be at list 2 input genomes"
        )
    )
    p.add_argument(
        "--hashsize",
        type=str,
        default="100M",
        help="Jellyfish initial hash size"
    )
    p.add_argument(
        "--in-memory",
        action="store_true",
        default=False,
        dest="in_memory",
        help="Keep Jellyfish counts in memory"
    )
    p.add_argument(
        "--k-min",
        type=int,
        default=4,
        dest="k_min",
        help="Minimum kmer size"
    )
    p.add_argument(
        "--k-max",
        type=int,
        dest="k_max",
        required=True,
        help="Maximum kmer size"
    )
    p.add_argument(
        "--lower",
        type=int,
        default=1,
        help="Do not let Jellyfish output kmers with count < --lower"
    )
    p.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Maximum number of CPUs to make it parallel"
    )
    p.add_argument(
        "--output",
        type=os.path.abspath,
        help="Path to the output file"
    )
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Maximum number of threads for Jellyfish"
    )
    return p.parse_args(args)


def count_kmers(
    genome: str,
    kmin: int = 4,
    kmax: int = 4,
    fast: bool = False,
    canonical: bool = False,
    lower: int = 1,
    hashsize: str = "100M",
    threads: int = 1
) -> Tuple[str, Dict[int, jf.Kmercount]]:
    """
    Run Jellyfish for counting kmers over an interval of kmer sizes on a input genome

    :param genome:      Path to the input genome
    :param kmin:        Minimum kmer size
    :param kmax:        Maximum kmer size
    :param fast:        Jellyfish one-pass calculation (faster)
    :param canonical:   Jellyfish count only canonical kmers
    :param lower:       Do not let Jellyfish output kmers with count < --lower
    :param hashsize:    Jellyfish initial hash size
    :param threads:     Make it parallel
    :return:            Dictionary with the list of kmers per kmer size
    """

    # Define a range of kmer sizes
    kmer_sizes = list(range(kmin, kmax + 1))

    datadict = dict()

    for kmer in kmer_sizes:
        datadict[kmer] = jf.Kmercount(
            genome,
            kmer,
            thread=threads,
            lower=lower,
            hashsize=hashsize,
            canonical=canonical,
            fast=fast
        )

    return genome, datadict


def filter(n: List[str], m: int) -> List[int]: 
    """
    Filter a list based on the occurrences of values

    :param n:   Input list
    :param m:   Filer factor
    :return:    Filtered list
    """

    c = collections.Counter(n)
    cr = set(n)
    result = list()

    for i in cr:
        if c[i] >= int(m * 0.5):
            result.append(i)

    return result


def par_acf(
    genome1: str,
    genome2: str,
    kmin: Optional[int] = None,
    kmax: Optional[int] = None,
    counts: Optional[Dict[str, Dict[int, jf.Kmercount]]] = None,
    fast: bool = False,
    canonical: bool = False,
    lower: int = 1,
    hashsize: str = "100M",
    threads: int = 1
) -> List[Tuple[int, int]]:
    """
    Make the computation of ACF parallel
    Run a instance of this function for each of the processes in a Pool

    :param genome1:     Genome #1
    :param genome2:     Genome #2
    :param kmin:        Minimum kmer size
    :param kmax:        Maximum kmer size
    :param counts:      Kmer counts for the input genomes
    :param fast:        Jellyfish one-pass calculation (faster)
    :param canonical:   Jellyfish count only canonical kmers
    :param lower:       Do not let Jellyfish output kmers with count < --lower
    :param hashsize:    Jellyfish initial hash size
    :param threads:     Maximum number of threads for Jellyfish
    :return:            ACF
    """

    res_acf = dict()

    for kmer in range(kmin, kmax + 1):
        if not counts:
            count_kmers_partial = partial(
                count_kmers,
                kmin=kmer,
                kmax=kmer,
                fast=fast,
                canonical=canonical,
                lower=lower,
                hashsize=hashsize,
                threads=threads
            )

            _, genome1_counts = count_kmers_partial(genome1)
            genome1_counts = genome1_counts[kmer]

            _, genome2_counts = count_kmers_partial(genome2)
            genome2_counts = genome2_counts[kmer]

            res_acf[kmer] = len(set(genome1_counts.keys()) & set(genome2_counts.keys()))

        else:
            res_acf[kmer] = len(set(counts[genome1][kmer].keys()) & set(counts[genome2][kmer].keys()))

    res_acf = sorted(res_acf.items(), key=itemgetter(0))

    return res_acf


def par_cre(
    genome: str,
    kmin: Optional[int] = None,
    kmax: Optional[int] = None,
    counts: Optional[Dict[str, Dict[int, jf.Kmercount]]] = None,
    fast: bool = False,
    canonical: bool = False,
    lower: int = 1,
    hashsize: str = "100M",
    threads: int = 1
) -> List[Tuple[int, float]]:
    """
    Make the computation of CRE parallel
    Run a instance of this function for each of the processes in a Pool

    :param genome:      Input genome
    :param kmin:        Minimum kmer size
    :param kmax:        Maximum kmer size
    :param counts:      Kmer counts for the input genomes
    :param fast:        Jellyfish one-pass calculation (faster)
    :param canonical:   Jellyfish count only canonical kmers
    :param lower:       Do not let Jellyfish output kmers with count < --lower
    :param hashsize:    Jellyfish initial hash size
    :param threads:     Maximum number of threads for Jellyfish
    :return:            CRE
    """

    res_cre = dict()

    def cal_re(a0: jf.Kmercount, a1: jf.Kmercount, a2: jf.Kmercount) -> float:
        """
        Calculate Relative Entropy (RE)

        :param a0:  Kmer count with kmer size
        :param a1:  Kmer count with kmer size - 1
        :param a2:  Kmer count with kmer size - 2
        :return:    Relative Entropy
        """

        result = 0.0
        rfactor = a0.sum
        lfactor = math.log((a1.sum ** 2) / (a2.sum * rfactor), 2)

        for key in a0.keys():
            realf = a0[key]
            left = a1[key[:-1]] if key[:-1] in a1 else 0
            right = a1[key[1:]] if key[1:] in a1 else 0
            below = a2[key[1:-1]] if key[1:-1] in a2 else 0

            if 0 not in (left, right, below):
                expectf = left * right / below
                result += max(0, realf / rfactor * (math.log(realf / expectf, 2) + lfactor))

        return result

    a0 = None
    a1 = None
    a2 = None

    count_kmers_partial = partial(
        count_kmers,
        fast=fast,
        canonical=canonical,
        lower=lower,
        hashsize=hashsize,
        threads=threads
    )

    for kmer in range(kmax, kmin - 1, -1):
        if not counts:
            if a0 is None and a1 is None and a2 is None:
                _, a0 = count_kmers_partial(genome, kmin=kmer, kmax=kmer)
                a0 = a0[kmer]

                _, a1 = count_kmers_partial(genome, kmin=(kmer - 1), kmax=(kmer - 1))
                a1 = a1[kmer - 1]

                _, a2 = count_kmers_partial(genome, kmin=(kmer - 2), kmax=(kmer - 2))
                a2 = a2[kmer - 2]

            else:
                a0 = a1
                a1 = a2

                _, a2 = count_kmers_partial(genome, kmin=(kmer - 2), kmax=(kmer - 2))
                a2 = a2[kmer - 2]
            
            res_cre[kmer] = cal_re(a0, a1, a2)

        else:
            res_cre[kmer] = cal_re(counts[genome][kmer], counts[genome][kmer - 1], counts[genome][kmer - 2])

        if kmer + 1 in res_cre:
            res_cre[kmer] += res_cre[kmer + 1]

    res_cre = sorted(res_cre.items(), key=itemgetter(0))

    return res_cre


def par_ofc(
    kmer: str,
    genomes: Optional[List[str]] = None,
    counts: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
    fast: bool = False,
    canonical: bool = False,
    lower: int = 1,
    hashsize: str = "100M",
    threads: int = 1
) -> Tuple[int, float]:
    """
    Make the computation of OCF parallel
    Run a instance of this function for each of the processes in a Pool

    :param kmer:        Kmer size
    :param genomes:     List of input genomes
    :param counts:      Kmer counts for the input genomes
    :param fast:        Jellyfish one-pass calculation (faster)
    :param canonical:   Jellyfish count only canonical kmers
    :param lower:       Do not let Jellyfish output kmers with count < --lower
    :param hashsize:    Jellyfish initial hash size
    :param threads:     Maximum number of threads for Jellyfish
    :return:            OCF
    """

    keys = list()

    count_kmers_partial = partial(
        count_kmers,
        kmin=kmer,
        kmax=kmer,
        fast=fast,
        canonical=canonical,
        lower=lower,
        hashsize=hashsize,
        threads=threads
    )

    for genome in genomes:
        if not counts:
            _, genome_counts = count_kmers_partial(genome)
            genome_counts = genome_counts[kmer]

        else:
            genome_counts = counts[genome][kmer]

        keys.extend(list(genome_counts.keys()))

    count_features = list(collections.Counter((collections.Counter(keys).values())).values())
    lnp = np.log2(count_features) - (kmer * 2)

    return kmer, np.sum(np.exp2(lnp) * lnp) * -1


def optimal_kmer_size(
    genomes: List[str],
    closely_related: bool = False,
    kmin: int = 4,
    kmax: Optional[int] = None,
    acf_cutoff: float = 0.1,
    cre_cutoff: float = 0.1,
    fast: bool = False,
    canonical: bool = False,
    lower: int = 1,
    hashsize: str = "100M",
    in_memory: bool = False,
    nproc: int = 1,
    threads: int = 1,
    output: Optional[str] = None
) -> None:
    """
    Compute the Cumulative Relative Entropy (CRE), Average Number of Common Features (ACF),
    and Observed Common Features on a seti of genomes over a kmer size interval.

    Eventually search for the optimal kmer size.

    :param genomes:         List with paths to the input genomes
    :param closely_related: Closely related genomes
    :param kmin:            Minimum kmer size (4 by default; cannot be lower than 4)
    :param kmax:            Maximum kmer size
    :param acf_cutoff:      Cutoff to use in selecting kmers whose ACF's are >= (cutoff * max(ACF))
    :param cre_cutoff:      Cutoff to use in selecting kmers whose CRE's are <= (cutoff * max(CRE))
    :param fast:            Jellyfish one-pass calculation (faster)
    :param canonical:       Jellyfish count only canonical kmers
    :param lower:           Do not let Jellyfish output kmers with count < --lower
    :param hashsize:        Jellyfish initial hash size
    :param in_memory:       Keep Jellyfish counts in memory
    :param nproc:           Maximum number of CPUs to make it parallel
    :param threads:         Maximum number of threads for Jellyfish
    :param output:          Path to the output file
    """

    if output:
        if os.path.isfile(output):
            raise Exception("Output file already exists")

        outfolder = os.path.dirname(output)

        if not os.path.isdir(outfolder):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), outfolder)

    if kmax <= kmin:
        raise ValueError("Invalid kmer size interval")

    counts = dict()

    if in_memory:
        print("Counting kmers in {} genomes".format(len(genomes)))

        if len(genomes) < 2:
            raise Exception("Not enough input genomes")

        with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(genomes)) as pbar:
            def progress(*args):
                pbar.update()

            count_kmers_partial = partial(
                count_kmers,
                kmin=(kmin - 2),
                kmax=kmax,
                fast=fast,
                canonical=canonical,
                lower=lower,
                hashsize=hashsize,
                threads=threads
            )

            jobs = [
                pool.apply_async(
                    count_kmers_partial,
                    args=(genome,),
                    callback=progress,
                )
                for genome in genomes
            ]

            for job in jobs:
                genome, genome_counts = job.get()
                counts[genome] = genome_counts

        if not counts:
            raise Exception("Something went wrong while counting kmer over the input genomes")

    if not closely_related:
        print("Computing Average number of Common Features (ACF) pairwise")

        acf_results = list()

        genomes_comb = combinations(genomes, 2)
        genomes_comb_len = int(math.factorial(len(genomes)) / (math.factorial(len(genomes) - 2) * 2))

        with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=genomes_comb_len) as pbar:
            def progress(*args):
                pbar.update()

            compute_acf_partial = partial(
                par_acf,
                kmin=kmin,
                kmax=kmax,
                counts=counts,
                fast=fast,
                canonical=canonical,
                lower=lower,
                hashsize=hashsize,
                threads=threads
            )

            jobs = [
                pool.apply_async(
                    compute_acf_partial,
                    args=(genome1, genome2,),
                    callback=progress,
                )
                for genome1, genome2 in genomes_comb
            ]

            for job in jobs:
                acf_results.append(job.get())

        acf_results = list(chain(*[result for result in acf_results]))

        acf_list = [x[1] for x in acf_results]
        maxacf = acf_cutoff * float(max(acf_list))
        acf_kmax = [x[0] for x in acf_results if x[1] >= maxacf]
        acf_kmax = filter(acf_kmax, genomes_comb_len)

        if acf_kmax:
            acf_kmax_max = max(acf_kmax)
            
            if acf_kmax_max != kmax:
                print("Redefining kmer size interval: --k-max={}".format(acf_kmax_max))

            kmax = acf_kmax_max

    else:
        acf_kmax = range(kmin, kmax + 1)

    print("Computing Cumulative Relative Entropy (CRE)")

    cre_results = list()

    with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(genomes)) as pbar:
        def progress(*args):
            pbar.update()
        
        compute_cre_partial = partial(
            par_cre,
            kmin=kmin,
            kmax=kmax,
            counts=counts,
            fast=fast,
            canonical=canonical,
            lower=lower,
            hashsize=hashsize,
            threads=threads
        )

        jobs = [
            pool.apply_async(
                compute_cre_partial,
                args=(genome,),
                callback=progress,
            )
            for genome in genomes
        ]

        for job in jobs:
            cre_results.append(job.get())

    cre_results = list(chain(*[result for result in cre_results]))
    
    cre_list = [x[1] for x in cre_results]
    maxval = 0 if not cre_list else cre_cutoff * float(max(cre_list))
    cre_kmin = [x[0] for x in cre_results if x[1] < maxval]
    cre_kmin = filter(cre_kmin, 0)

    if cre_kmin:
        cre_kmin_min = min(cre_kmin)
            
        if cre_kmin_min != kmin:
            print("Redefining kmer size interval: --k-min={}".format(cre_kmin_min))

        kmin = cre_kmin_min

    print("Computing Observed Common Features (OCF)")

    ofc_results = dict()

    kmer_sizes = list(range(kmin, kmax + 1))

    with mp.Pool(processes=nproc) as pool, tqdm.tqdm(total=len(kmer_sizes)) as pbar:
        def progress(*args):
            pbar.update()

        compute_ofc_partial = partial(
            par_ofc,
            genomes=genomes,
            counts=counts,
            fast=fast,
            canonical=canonical,
            lower=lower,
            hashsize=hashsize,
            threads=threads
        )

        jobs = [
            pool.apply_async(
                compute_ofc_partial,
                args=(kmer,),
                callback=progress,
            )
            for kmer in kmer_sizes
        ]

        for job in jobs:
            kmer, res = job.get()
            ofc_results[kmer] = res

    ofc_results = sorted(ofc_results.items(), key=itemgetter(0))

    ofc_list = [x[1] for x in ofc_results]

    if len(ofc_list) != 0:
        ofc_max = max(ofc_list)

    else:
        ofc_max = 1e7

    k_max = [x for x, y in ofc_results if y == ofc_max]
    ofc_possible_kmers = [x[0] for x in ofc_results if x[0] >= k_max[0]]

    opt = list(set(cre_kmin) & set(acf_kmax) & set(ofc_possible_kmers))

    if opt:
        out_message = "The recommended choice of kmer is: {}".format(min(opt))

    else:
        out_message = ("The recommended choice of kmer does not lie within the given interval "
                       "from 4 to the largest kmer length you selected")

    print(out_message, file=open(output, "w+") if output else None)


def run(args) -> None:
    # Load command line parameters
    args = read_params(args)

    t0 = time.time()

    optimal_kmer_size(
        [line.strip() for line in open(args.filenames).readlines() if line.strip()],
        closely_related=args.closely_related,
        kmin=args.k_min,
        kmax=args.k_max,
        acf_cutoff=args.acf_cutoff,
        cre_cutoff=args.cre_cutoff,
        fast=args.fast,
        canonical=args.canonical,
        lower=args.lower,
        hashsize=args.hashsize,
        in_memory=args.in_memory,
        nproc=args.nproc,
        threads=args.threads,
        output=args.output
    )

    t1 = time.time()

    print("Total elapsed time {}s".format(int(t1 - t0)))


if __name__ == "__main__":
    run(sys.argv)
