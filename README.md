# KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics

<img src="https://github.com/natapol/kitsune/blob/master/logoKITSUNE.png" width="40%">

[![PyPI version](https://badge.fury.io/py/kitsune.svg)](https://badge.fury.io/py/kitsune)
![Upload Python Package](https://github.com/natapol/kitsune/workflows/Upload%20Python%20Package/badge.svg)

KITSUNE is a toolkit for evaluation of the length of k-mer in a given genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many applications including biological sequence comparison. However, selection of an appropriate k-mer length to obtain a good information content for comparison is normally overlooked. The optimum k-mer length is a prerequsite to obtain biological meaningful genomic distance for assesment of phylogenetic relationships. Therefore, we have developed KITSUNE to aid k-mer length selection process in a systematic way, based on a three-steps aproach described in [Viral Phylogenomics Using an Alignment-Free Method: A Three-Step Approach to Determine Optimal Length of k-mer](https://doi.org/10.1038/srep40712).

KITSUNE will calculte the three matrices across considered k-mer range:

1. Cumulative Relative Entropy (CRE)
1. Average number of Common Features (ACF)
1. Observed Common Features (OCF)

Moreover, KITSUNE also provides various genomic distance calculations from the k-mer frequency vectors that can be used for species identification or phylogenomic tree construction.  

> **Note**:
> If you use KITSUNE in your research, please cite:
> [KITSUNE: A Tool for Identifying Empirically Optimal K-mer Length for Alignment-free Phylogenomic Analysis](https://doi.org/10.3389/fbioe.2020.556413)


## Installation

Kitsune is developed under python version 3 environment. We recommend users use python >= v3.5.

Requirement packages: scipy >= 0.18.1, numpy >= 1.1.0, tqdm >= 4.32

Kitsune also requires [Jellyfish](https://doi.org/10.1093/bioinformatics/btr011) for k-mer counting as external software dependency. Thus, you need to install it before running the tool: [https://github.com/gmarcais/Jellyfish](https://github.com/gmarcais/Jellyfish)

### Install with pip

```bash
pip install kitsune
```

### Install from source

```bash
# Clone the GitHub repository
git clone https://github.com/natapol/kitsune

# Move to the kitsune folder
cd kitsune/

# Install
python setup.py install
```

## Usage

### Overview of kitsune

command for listing help

```bash
$ kitsune --help

usage: kitsune <command> [<args>]

Available commands:
  acf      Compute average number of common features between signatures
  cre      Compute cumulative relative entropy
  dmatrix  Compute distance matrix
  kopt     Compute recommended choice (optimal) of kmer within a given kmer interval for a set of genomes using the cre, acf and ofc
  ofc      Compute observed feature frequencies

Use --help in conjunction with one of the commands above for a list of available options (e.g. kitsune acf --help)
```

### Calculate CRE, ACF, and OFC value for specific kmer

Kitsune provides three commands to calculate an appropiate k-mer using CRE, ACF, and OCF:

#### Calculate CRE

```bash
$ kitsune cre -h

usage: kitsune (cre) [-h] --filename FILENAME [--fast] [--canonical] -ke KEND
                     [-kf KFROM] [-t THREAD] [-o OUTPUT]

Calculate k-mer from cumulative relative entropy of all genomes

optional arguments:
  -h, --help            show this help message and exit
  --filename FILENAME   A genome file in fasta format (default: None)
  --fast                Jellyfish one-pass calculation (faster) (default:
                        False)
  --canonical           Jellyfish count only canonical mer (default: False)
  -ke KEND, --kend KEND
                        Last k-mer (default: None)
  -kf KFROM, --kfrom KFROM
                        Calculate from k-mer (default: 4)
  -t THREAD, --thread THREAD
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
```

#### Calculate ACF

```bash
$ kitsune acf -h

usage: kitsune (acf) [-h] --filenames FILENAMES [FILENAMES ...] [--fast]
                     [--canonical] -k KMERS [KMERS ...] [-t THREAD]
                     [-o OUTPUT]

Calculate an average number of common features pairwise between one genome
against others

optional arguments:
  -h, --help            show this help message and exit
  --filenames FILENAMES [FILENAMES ...]
                        Genome files in fasta format (default: None)
  --fast                Jellyfish one-pass calculation (faster) (default:
                        False)
  --canonical           Jellyfish count only canonical mer (default: False)
  -k KMERS [KMERS ...], --kmers KMERS [KMERS ...]
                        Have to state before (default: None)
  -t THREAD, --thread THREAD
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
```

#### Calculate OFC

```bash
$ kitsune ofc -h

usage: kitsune (ofc) [-h] --filenames FILENAMES [FILENAMES ...] [--fast]
                     [--canonical] -k KMERS [KMERS ...] [-t THREAD]
                     [-o OUTPUT]

Calculate an observe feature frequency

optional arguments:
  -h, --help            show this help message and exit
  --filenames FILENAMES [FILENAMES ...]
                        Genome files in fasta format (default: None)
  --fast                Jellyfish one-pass calculation (faster) (default:
                        False)
  --canonical           Jellyfish count only canonical mer (default: False)
  -k KMERS [KMERS ...], --kmers KMERS [KMERS ...]
  -t THREAD, --thread THREAD
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
```

#### General Example

```bash
kitsune cre --filename genome1.fna -kf 5 -ke 10
kitsune acf --filenames genome1.fna genome2.fna -k 5
kitsune ofc --filenames genome_fasta/* -k 5
```

### Calculate genomic distance at specific k-mer from kmer frequency vectors of two of genomes

Kitsune provides a commands to calculate genomic distance using different distance estimation method. Users can assess the impact of a selected k-mer length on the genomic distnace of choice below.

distance option        | name
---------------------- | ----
braycurtis             | Bray-Curtis distance
canberra               | Canberra distance
chebyshev              | Chebyshev distance
cityblock              | City Block (Manhattan) distance
correlation            | Correlation distance
cosine                 | Cosine distance
euclidean              | Euclidean distance
jensenshannon          | Jensen-Shannon distance
sqeuclidean            | Squared Euclidean distance
dice                   | Dice dissimilarity
hamming                | Hamming distance
jaccard                | Jaccard-Needham dissimilarity
kulsinski              | Kulsinski dissimilarity
rogerstanimoto         | Rogers-Tanimoto dissimilarity
russellrao             | Russell-Rao dissimilarity
sokalmichener          | Sokal-Michener dissimilarity
sokalsneath            | Sokal-Sneath dissimilarity
yule                   | Yule dissimilarity
mash                   | MASH distance
jsmash                 | MASH Jensen-Shannon distance
jaccarddistp           | Jaccard-Needham dissimilarity Probability
euclidean_of_frequency | Euclidean distance of Frequency


Kitsune provides a choice of distance transformation proposed by [Fan et.al](https://doi.org/10.1186/s12864-015-1647-5).

### Calculate a distance matrix

```bash
$ kitsune dmatrix -h

usage: kitsune (dmatrix) [-h] [--filenames [FILENAMES [FILENAMES ...]]]
                         [--fast] [--canonical] -k KMER [-i INPUT] [-o OUTPUT]
                         [-t THREAD] [--transformed]
                         [-d {braycurtis,canberra,jsmash,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,kulsinsk,matching,rogerstanimoto,russellrao,sokalmichener,sokalsneath,sqeuclidean,yule,mash,jaccarddistp}]
                         [-f FORMAT]

Calculate a distance matrix

optional arguments:
  -h, --help            show this help message and exit
  --filenames [FILENAMES [FILENAMES ...]]
                        Genome files in fasta format (default: None)
  --fast                Jellyfish one-pass calculation (faster) (default:
                        False)
  --canonical           Jellyfish count only canonical mer (default: False)
  -k KMER, --kmer KMER
  -i INPUT, --input INPUT
                        List of genome files in txt (default: None)
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
  -t THREAD, --thread THREAD
  --transformed
  -d {braycurtis,canberra,jsmash,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,kulsinsk,matching,rogerstanimoto,russellrao,sokalmichener,sokalsneath,sqeuclidean,yule,mash,jaccarddistp}, --distance {braycurtis,canberra,jsmash,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,kulsinsk,matching,rogerstanimoto,russellrao,sokalmichener,sokalsneath,sqeuclidean,yule,mash,jaccarddistp}
  -f FORMAT, --format FORMAT
```

Example of choosing distance option:

```bash
kitsune dmatrix --filenames genome1.fna genome2.fna -k 11 -d jaccard --canonical --fast -o output.txt
kitsune dmatrix --filenames genome1.fna genome2.fna -k 11 -d hensenshannon --canonical --fast -o output.txt
```

### Find optimum k-mer from a given set of genomes

Kitsune provides a wrap-up comand to find optimum k-mer length for a given set of genome within a given kmer interval.

```bash
$ kitsune kopt -h

usage: kitsune (kopt) [-h] [--acf-cutoff ACF_CUTOFF] [--canonical]
                      [--closely-related] [--cre-cutoff CRE_CUTOFF] [--fast]
                      --filenames FILENAMES [--hashsize HASHSIZE]
                      [--in-memory] [--k-min K_MIN] --k-max K_MAX
                      [--lower LOWER] [--nproc NPROC] [--output OUTPUT]
                      [--threads THREADS]

Optimal kmer size selection for a set of genomes using Average number of
Common Features (ACF), Cumulative Relative Entropy (CRE), and Observed Common
Features (OCF). Example: kitsune kopt --filenames genomeList.txt --k-min 4
--k-max 12 --canonical --fast

optional arguments:
  -h, --help            show this help message and exit
  --acf-cutoff ACF_CUTOFF
                        Cutoff to use in selecting kmers whose ACFs are >=
                        (cutoff * max(ACF)) (default: 0.1)
  --canonical           Jellyfish count only canonical kmers (default: False)
  --closely-related     Use in case of closely related genomes (default:
                        False)
  --cre-cutoff CRE_CUTOFF
                        Cutoff to use in selecting kmers whose CREs are <=
                        (cutoff * max(CRE)) (default: 0.1)
  --fast                Jellyfish one-pass calculation (faster) (default:
                        False)
  --filenames FILENAMES
                        Path to the file with the list of genome files paths.
                        There should be at list 2 input genomes (default:
                        None)
  --hashsize HASHSIZE   Jellyfish initial hash size (default: 100M)
  --in-memory           Keep Jellyfish counts in memory (default: False)
  --k-min K_MIN         Minimum kmer size (default: 4)
  --k-max K_MAX         Maximum kmer size (default: None)
  --lower LOWER         Do not let Jellyfish output kmers with count < --lower
                        (default: 1)
  --nproc NPROC         Maximum number of CPUs to make it parallel (default:
                        1)
  --output OUTPUT       Path to the output file (default: None)
  --threads THREADS     Maximum number of threads for Jellyfish (default: 1)
```

### Example dataset

First download the example files. [Download](examples/example_viral_genomes.zip)

```bash
kitsune kopt --filenames genome_list --k-min 6 --k-max 21 --canonical --fast --threads 4 --nproc 2 --output out.txt
```

> :warning: _Please be aware that this command will use big computational resources when large number of genomes and/or large genome size are used as the input._
