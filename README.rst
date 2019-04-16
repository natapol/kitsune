kitsune
--------

# Kitsune
Kitsune is a toolkit for evaluate the lenght of optimal k-mer in dataset.

## Basics
Kitsune provides many methods to evaluate the amount of information by given k-mer.
This quantitative measurement could be use to calculate the amount of need k-mer

## Installation
Clone the repository and install it throught pip
```
pip install kitsune
```

## Usage

### Calculate k-mer from CRE, ACF, and OFC value
Kitsune provides three commands for calculate an appropiate k-mer using CRE, ACF, and OCF.

```
kitsune cre genome_fasta/* -ks 5 -ke 10\n
kitsune acf genome_fasta/* -ks 5 -ke 10\n
kitsune ocf genome_fasta/* -ks 5 -ke 10
```
