# Ksiga
KSiga is a toolkit for evaluate the lenght of optimal k-mer in dataset.

## Basics
KSiga provides many methods to evaluate the amount of information by given k-mer.
This quantitative measurement could be use to calculate the amount of need k-mer

## Installation
Clone the repository and install it throught pip
```
git clone https://github.com/yumyai/ksiga.git
cd ksiga
pip install -e .
```

## Usage

### Calculate k-mer from CRE, ACF, and OFC value
KSiga provides three commands for calculate an appropiate k-mer using CRE, ACF, and OCF.

```
ksiga cre genome_fasta/* -ks 5 -ke 10
ksiga acf genome_fasta/* -ks 5 -ke 10
ksiga ocf genome_fasta/* -ks 5 -ke 10
```

