KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics
=======================================================================

.. figure:: https://github.com/natapol/kitsune/blob/master/logoKITSUNE.png?v&s=200
   :alt: KISUNE

   KISUNE

KITSUNE is a toolkit for evaluation of the lenght of k-mer in a given
genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many
applications including biological sequence comparison. However,
selection of an appropriate k-mer length to obtain a good information
content for comparison is normally overlooked. Therefore, we have
developed KITSUNE to aid k-mer length selection process based on a three
steps aproach described in `Viral Phylogenomics Using an Alignment-Free
Method: A Three-Step Approach to Determine Optimal Length of
k-mer <https://www.nature.com/articles/srep40712>`__.

KITSUNE uses Jellyfish software
`Jellyfish <https://academic.oup.com/bioinformatics/article/27/6/764/234905>`__
for k-mer counting. Thanks to Jellyfish developer.

KITSUNE will calculte the three matrices across considered k-emer range
:

1. Cumulative Relative Entropy (CRE)
2. Averrage number of Common Feature (ACF)
3. Obserbed Common Feature (OCF)

Moreverver, KITSUNE also provides various genomic distance calculations
from the k-mer frequnce vectors that can be used for species
identifiction or phylogenomic tree construction.

If you use KITSUNE in your research, please cite:
`Reference <https://github.com/natapol/kitsune>`__

Installation
------------

Clone the repository and install it throught pip

.. code:: bash

   pip install kitsune

Usage
-----

Calculate CRE, ACF, and OFC value for specific kmer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides three commands to calculate an appropiate k-mer using
CRE, ACF, and OCF.

.. code:: bash

   kitsune cre genome_fasta/* -ks 5 -ke 10
   kitsune acf genome_fasta/* -ks 5 -ke 10
   kitsune ocf genome_fasta/* -ks 5 -ke 10

Calculate genomic distance at specific k-mer from kmer frequency vectors of two of genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides a commands to calculate genomic distance using
different distance estimation method.

=============== =========================================
distance option name
=============== =========================================
braycurtis      Bray-Curtis distance
canberra        Canberra distance
chebyshev       Chebyshev distance
cityblock       City Block (Manhattan) distance
correlation     Correlation distance
cosine          Cosine distance
euclidean       Euclidean distance
jensenshannon   Jensen-Shannon distance
sqeuclidean     Squared Euclidean distance
dice            Dice dissimilarity
hamming         Hamming distance
jaccard         Jaccard-Needham dissimilarity
kulsinski       Kulsinski dissimilarity
rogerstanimoto  Rogers-Tanimoto dissimilarity
russellrao      Russell-Rao dissimilarity
sokalmichener   Sokal-Michener dissimilarity
sokalsneath     Sokal-Sneath dissimilarity
yule            Yule dissimilarity
mash            MASH distance
jsmash          MASH Jensen-Shannon distance
jaccarddistp    Jaccard-Needham dissimilarity Probability
=============== =========================================

.. code:: bash

   kitsune dmatrix genome1.fna genome2.fna -k 17 -d jaccard --canonical --fast -o output.txt
   kitsune dmatrix genome1.fna genome2.fna -k 17 -d hensenshannon --canonical --fast -o output.txt

Find optimum k-mer from a given set of genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides a comand to find optimum k-mer length in agiven set of
genome.

First download the example
files.\ `Download <%22https://github.com/natapol/kitsune/blob/master/examaple_viral_genomes.zip%22>`__

Then use kitsune kopt command

-i : path to list of genome files

-ks: The smallest kmer-length to consider

-kl: The largest kmer-length to consider

-o: output file

\**Please be aware that this comand will use big computational resources
when large number of genomes and/or large genome size are used as the
input.

.. code:: bash

   kitsune kopt -i genome_list -ks 7 -kl 15 --canonical --fast -o output.txt
