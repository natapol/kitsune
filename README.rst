KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics
=======================================================================

|PyPI version| |Upload Python Package|

KITSUNE is a toolkit for evaluation of the length of k-mer in a given
genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many
applications including biological sequence comparison. However,
selection of an appropriate k-mer length to obtain a good information
content for comparison is normally overlooked. The optimum k-mer length
is a prerequsite to obtain biological meaningful genomic distance for
assesment of phylogenetic relationships. Therefore, we have developed
KITSUNE to aid k-mer length selection process in a systematic way, based
on a three-steps aproach described in `Viral Phylogenomics Using an
Alignment-Free Method: A Three-Step Approach to Determine Optimal Length
of k-mer <https://www.nature.com/articles/srep40712>`__.

KITSUNE uses Jellyfish software for k-mer counting. Thanks to Jellyfish
developer.
`Citation <https://academic.oup.com/bioinformatics/article/27/6/764/234905>`__

KITSUNE will calculte the three matrices across considered k-mer range:

1. Cumulative Relative Entropy (CRE)
2. Averrage number of Common Feature (ACF)
3. Obserbed Common Feature (OCF)

Moreover, KITSUNE also provides various genomic distance calculations
from the k-mer frequency vectors that can be used for species
identification or phylogenomic tree construction.

If you use KITSUNE in your research, please cite: KITSUNE: A Tool for
Identifying Optimal K-mer Length for Alignment-free Phylogenomic
Analysis `Reference <https://github.com/natapol/kitsune>`__

Installation
------------

Kitsune is developed under python version 3 environment. We recommend
users use python >= v3.5.

Requirement packages:

biopython >= 1.68, scipy >= 0.18.1, numpy >= 1.1.0, tqdm >= 4.32

pip
~~~

.. code:: bash

   pip install kitsune

Clone from github
~~~~~~~~~~~~~~~~~

.. code:: bash

   git clone https://github.com/natapol/kitsune
   cd kitsune/
   python nstall setup.py

Usage
-----

Overview of kitsune
~~~~~~~~~~~~~~~~~~~

command for listing help

.. code:: bash

   $ kitsune --help

   usage: kitsune <command> [<args>]

   Commands can be:
   cre <filename>                    Compute cumulative relative entropy.
   acf <filenames>                   Compute average number of common feature between signatures.
   ofc <filenames>                   Compute observed feature frequencies.
   kopt <filenames>                  Compute recommended choice (optimal) of kmer within a given kmer interval for a set of genomes using the cre, acf and ofc.
   dmatrix <filenames>               Compute distance matrix.

Calculate CRE, ACF, and OFC value for specific kmer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides three commands to calculate an appropiate k-mer using
CRE, ACF, and OCF:

Calculate CRE
~~~~~~~~~~~~~

.. code:: bash

   $ kitsune cre -h
   usage: kitsune [-h] [--fast] [--canonical] -ke KEND [-kf KFROM] [-t THREAD]
                  [-o OUTPUT]
                  filename

   Calculate k-mer from cumulative relative entropy of all genomes

   positional arguments:
     filename              a genome file in fasta format

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -ke KEND, --kend KEND
                           last k-mer
     -kf KFROM, --kfrom KFROM
                           Calculate from k-mer
     -t THREAD, --thread THREAD
     -o OUTPUT, --output OUTPUT
                           output filename

Calculate ACF
~~~~~~~~~~~~~

.. code:: bash

   $ kitsune acf -h
   usage: kitsune [-h] [--fast] [--canonical] -k KMERS [KMERS ...] [-t THREAD]
                  [-o OUTPUT]
                  filenames [filenames ...]

   Calculate average number of common feature

   positional arguments:
     filenames             genome files in fasta format

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -k KMERS [KMERS ...], --kmers KMERS [KMERS ...]
                           have to state before
     -t THREAD, --thread THREAD
     -o OUTPUT, --output OUTPUT
                           output filename

Calculate OFC
~~~~~~~~~~~~~

.. code:: bash

   $ kitsune ofc -h
   usage: kitsune [-h] [--fast] [--canonical] -k KMERS [KMERS ...] [-t THREAD]
                  [-o OUTPUT]
                  filenames [filenames ...]

   Calculate observe feature occurrence

   positional arguments:
     filenames             genome files in fasta format

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -k KMERS [KMERS ...], --kmers KMERS [KMERS ...]
     -t THREAD, --thread THREAD
     -o OUTPUT, --output OUTPUT
                           output filename

General Example
~~~~~~~~~~~~~~~

.. code:: bash

   kitsune cre genome1.fna -kf 5 -ke 10
   kitsune acf genome1.fna genome2.fna -k 5
   kitsune ofc genome_fasta/* -k 5

Calculate genomic distance at specific k-mer from kmer frequency vectors of two of genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides a commands to calculate genomic distance using
different distance estimation method. Users can assess the impact of a
selected k-mer length on the genomic distnace of choice below.

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

Kitsune provides a choice of distance transformation proposed by `Fan
et.al <https://doi.org/10.1186/s12864-015-1647-5>`__.

Calculate a distance matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   $ kitsune dmatrix -h
   usage: kitsune [-h] [--fast] [--canonical] -k KMER [-i INPUT] [-o OUTPUT]
                  [-t THREAD] [--transformed] [-d DISTANCE] [-f FORMAT]
                  [filenames [filenames ...]]

   Calculate a distance matrix

   positional arguments:
     filenames             genome files in fasta format

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -k KMER, --kmer KMER
     -i INPUT, --input INPUT
                           list of genome files in txt
     -o OUTPUT, --output OUTPUT
                           output filename
     -t THREAD, --thread THREAD
     --transformed
     -d DISTANCE, --distance DISTANCE
                           braycurtis, canberra, jsmash, chebyshev, cityblock,
                           correlation, cosine (default), dice, euclidean,
                           hamming, jaccard, kulsinsk, matching, rogerstanimoto,
                           russellrao, sokalmichener, sokalsneath, sqeuclidean,
                           yule, mash, jaccarddistp
     -f FORMAT, --format FORMAT

Example of choosing distance option:

.. code:: bash

   kitsune dmatrix genome1.fna genome2.fna -k 11 -d jaccard --canonical --fast -o output.txt
   kitsune dmatrix genome1.fna genome2.fna -k 11 -d hensenshannon --canonical --fast -o output.txt

Find optimum k-mer from a given set of genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides a wrap-up comand to find optimum k-mer length for a
given set of genome within a given kmer interval.

.. code:: bash

   $ kitsune kopt -h
   usage: kitsune [-h] [--fast] [--canonical] -kl KLARGE [-o OUTPUT]
                  [--closely_related] [-x CRE_CUTOFF] [-y ACF_CUTOFF] [-t THREAD]
                  filenames

   Example: kitsune kopt genomeList.txt -kl 15 --canonical --fast -t 4 -o out.txt

   positional arguments:
     filenames             A file that list the path to all genomes(fasta format)
                           with extension as (.txt,.csv,.tab) or no extension

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -kl KLARGE, --klarge KLARGE
                           largest k-mer length to consider, note: the smallest
                           kmer length is 4
     -o OUTPUT, --output OUTPUT
                           output filename
     --closely_related     For closely related set of genomes, use this option
     -x CRE_CUTOFF, --cre_cutoff CRE_CUTOFF
                           cutoff to use in selecting kmers whose cre's are <=
                           (cutoff * max(cre)), Default = 10 percent, ie x=0.1
     -y ACF_CUTOFF, --acf_cutoff ACF_CUTOFF
                           cutoff to use in selecting kmers whose acf's are >=
                           (cutoff * max(acf)), Default = 10 percent, ie y=0.1
     -t THREAD, --thread THREAD
                           Number of threads (integer)

Example dataset
~~~~~~~~~~~~~~~

First download the example files.
`Download <https://github.com/natapol/kitsune/blob/master/examaple_viral_genomes.zip>`__

.. code:: bash

   kitsune kopt genomeList.txt -kl 15 --canonical --fast -t 4 -o out.txt

\**Please be aware that this command will use big computational
resources when large number of genomes and/or large genome size are used
as the input.

.. |PyPI version| image:: https://badge.fury.io/py/kitsune.svg
   :target: https://badge.fury.io/py/kitsune
.. |Upload Python Package| image:: https://github.com/natapol/kitsune/workflows/Upload%20Python%20Package/badge.svg
