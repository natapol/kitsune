add logo here

# KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics

KITSUNE is a toolkit for evaluate the lenght of optimal k-mer length in a given genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many applications including biological sequence comparison. However, selection of an appropriate k-mer length to obtain a good information is normally overlooked. Therefore, we have develop KITSUNE to aid k-mer length selection process base on a three steps aproach (doi: 10.1038/srep40712 (2017). 

KITSUNE will calculte three matrices across considered k-emer range :
(1)	Cumulative Relative Entropy (CRE)
(2) Averrage number of Common Feature (ACF)
(3) Obserb Feature Occurance (OCF) 

Mooverver, KITSUNE also provide various genomic distance calculations from the k-mer frequnce vector that can be used for species identifiction or phylogenomic tree construction  

## Installation
Clone the repository and install it throught pip
```
pip install kitsune .
```

## Usage

### Calculate k-mer from CRE, ACF, and OFC value
Kitsune provides three commands for calculate an appropiate k-mer using CRE, ACF, and OCF.

```
kitsune cre genome_fasta/* -ks 5 -ke 10
kitsune acf genome_fasta/* -ks 5 -ke 10
kitsune ocf genome_fasta/* -ks 5 -ke 10

calculate distance .....
```

### Example

xyz
