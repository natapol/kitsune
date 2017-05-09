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
### Create k-mer profiles.
To create a k-mer profiles.
```
wget https://github.com/yumyai/ksiga/blob/devel/examples.tar.gz?raw=true  # Download example dataset
tar -xzf examples.tar.gz
# This will create an index in folder ex-index
ksiga index examples -o ex-index
```

### Calculate k-mer from CRE, ACF, and OFC value
KSiga provides three commands for calculate an appropiate k-mer using CRE, ACF, and OCF.

```
ksiga cre_kmer ex-index/* -ks 5 -ke 10
ksiga acf_kmer ex-index/* -ks 5 -ke 10
ksiga ocf_kmer ex-index/* -ks 5 -ke 10
```

Each of these command will print out the optimal k-mer.

### Standalone CRE calculation.
Calculating CRE are very computing intensive. We planned to provide a way to run in on cluster soon.
