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
wget examples.tar.gz # Download example dataset
tar -xzf examples.tar.gz
# This will create an index in folder ex-index
ksiga index examples -o ex-index
```

### Calculate CRE, ACF, and OFC value
KSiga provides three commands for calculate OF

```
ksiga rel
ksiga cre
ksiga acf
ksiga ocf
```
