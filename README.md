![KISUNE](https://github.com/natapol/kitsune/blob/master/logoKITSUNE.png?v&s=200)

# KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics

KITSUNE is a toolkit for evaluation of the lenght of k-mer in a given genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many applications including biological sequence comparison. However, selection of an appropriate k-mer length to obtain a good information content for comparison is normally overlooked. Therefore, we have developed KITSUNE to aid k-mer length selection process based on a three steps aproach (https://www.nature.com/articles/srep40712). 

KITSUNE uses Jellyfish software (doi:10.1093/bioinformatics/btr011) for k-mer counting. Thanks to Jellyfish developer.

KITSUNE will calculte the three matrices across considered k-emer range :
(1)	Cumulative Relative Entropy (CRE),
(2) Averrage number of Common Feature (ACF),
(3) Obserbed Common Feature (OCF) 

Moreverver, KITSUNE also provides various genomic distance calculations from the k-mer frequnce vectors that can be used for species identifiction or phylogenomic tree construction.  

If you use KITSUNE in your research, please cite:
xyx

## Installation
Clone the repository and install it throught pip
```
pip install kitsune .
```

## Usage

### Calculate CRE, ACF, and OFC value for specific kmer.
Kitsune provides three commands to calculate an appropiate k-mer using CRE, ACF, and OCF.

```
kitsune cre genome_fasta/* -ks 5 -ke 10
kitsune acf genome_fasta/* -ks 5 -ke 10
kitsune ocf genome_fasta/* -ks 5 -ke 10
```

### Calculate genomic distance at specific k-mer from kmer frequency vectors of two of genomes.

Kitsune provides a commands to calculate genomic distance using different distance estimation method.

```
kitsune dmatrix genome1.fna genome2.fna -k 17 -d jaccard --canonical --fast -o output.txt
kitsune dmatrix genome1.fna genome2.fna -k 17 -d hensenshannon --canonical --fast -o output.txt
```

### Find optimum k-mer from a given set of genome.

Kitsune provides a comand to find optimum k-mer length in agiven set of genome. 

First download the example files.[Download]("https://github.com/natapol/kitsune/blob/master/examaple_viral_genomes.zip") 

Then use kitsune kopt command

-i : path to list of genome files

-ks: The smallest kmer-length to consider

-kl: The largest kmer-length to consider

-o: output file

**Please be aware that this comand will use big computational resources when large number of genomes and/or lrage geenome siz are used as the input.    

```
kitsune kopt -i genome_list -ks 7 -kl 15 --canonical --fast -o output.txt
```
