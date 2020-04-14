![KISUNE](https://github.com/natapol/kitsune/blob/master/logoKITSUNE.png){:height="50%" width="50%"}

# KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics

KITSUNE is a toolkit for evaluation of the lenght of k-mer in a given genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many applications including biological sequence comparison. However, selection of an appropriate k-mer length to obtain a good information content for comparison is normally overlooked. Therefore, we have developed KITSUNE to aid k-mer length selection process based on a three steps aproach (https://www.nature.com/articles/srep40712). 

KITSUNE uses Jellyfish software (doi:10.1093/bioinformatics/btr011) for k-mer counting. Thanks to Jellyfish developer.

KITSUNE will calculte the three matrices across considered k-emer range :
(1)	Cumulative Relative Entropy (CRE)
(2) Averrage number of Common Feature (ACF)
(3) Obserb Feature Occurance (OCF) 

Moreverver, KITSUNE also provides various genomic distance calculations from the k-mer frequnce vectors that can be used for species identifiction or phylogenomic tree construction.  

If you use KITSUNE in your research, please cite:
xyx

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
