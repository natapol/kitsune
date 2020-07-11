KITSUNE: K-mer-length Iterative Selection for UNbiased Ecophylogenomics
=======================================================================

|PyPI version| |Upload Python Package|

KITSUNE is a toolkit for evaluation of the lenght of k-mer in a given
genome dataset for alignment-free phylogenimic analysis.

K-mer based approach is simple and fast yet has been widely used in many
applications including biological sequence comparison. However,
selection of an appropriate k-mer length to obtain a good information
content for comparison is normally overlooked. The optimum k-mer length
is a prerequsite to obtain biological menaingful genomic distance for
assesment of phylogenetic relationships. Therefore, we have developed
KITSUNE to aid k-mer length selection process in a systematic way, based
on a three-steps aproach described in `Viral Phylogenomics Using an
Alignment-Free Method: A Three-Step Approach to Determine Optimal Length
of k-mer <https://www.nature.com/articles/srep40712>`__.

KITSUNE uses Jellyfish software for k-mer counting. Thanks to Jellyfish
developer.
`Citation <https://academic.oup.com/bioinformatics/article/27/6/764/234905>`__

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

Kitsune is developed under python version 3 environment. We recomend
python v3.8……. Requirement packages scipy>=0.18.1, numpy>=1.1.0,
qdm>=4.32

.. code:: bash

   pip install kitsune

Usage
-----

Overview of kitsune
-------------------

.. code:: bash

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
                           

Example
~~~~~~~

.. code:: bash

   kitsune cre genome_fasta/* -kf 5 -ke 10
   kitsune acf genome_fasta/* -k 5
   kitsune ofc genome_fasta/* -k 5 

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

Calculate a distance matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   positional arguments:
     filenames             genome files in fasta format

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer (use for raw read
                           count)
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

   kitsune dmatrix genome1.fna genome2.fna -k 17 -d jaccard --canonical --fast -o output.txt
   kitsune dmatrix genome1.fna genome2.fna -k 17 -d hensenshannon --canonical --fast -o output.txt

Find optimum k-mer from a given set of genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kitsune provides a comand to find optimum k-mer length for a given set
of genome within a given kmer interval.

.. code:: bash

   usage: kitsune [-h] [--fast] [--canonical] -kl KLARGE [-o OUTPUT] [--closely_related] [-x CRE_CUTOFF] [-y ACF_CUTOFF] [-t THREAD] filenames

   Example: kitsune kopt genomeList.txt -kl 15 --canonical --fast -t 4 -o out.txt

   positional arguments:
     filenames             A file that list the path to all genomes(fasta format) with extension as (.txt,.csv,.tab) or no extension

   optional arguments:
     -h, --help            show this help message and exit
     --fast                Jellyfish one-pass calculation (faster)
     --canonical           Jellyfish count only canonical mer
     -kl KLARGE, --klarge KLARGE
                           largest k-mer length to consider, note: the smallest kmer length is 4
     -o OUTPUT, --output OUTPUT
                           output filename
     --closely_related     For closely related set of genomes, use this option
     -x CRE_CUTOFF, --cre_cutoff CRE_CUTOFF
                           cutoff to use in selecting kmers whose cre's are <= (cutoff * max(cre)), Default = 10 percent, ie x=0.1
     -y ACF_CUTOFF, --acf_cutoff ACF_CUTOFF
                           cutoff to use in selecting kmers whose acf's are >= (cutoff * max(acf)), Default = 10 percent, ie y=0.1
     -t THREAD, --thread THREAD
                           Number of threads (integer)

First download the example
files.\ `Download <%22https://github.com/natapol/kitsune/blob/master/examaple_viral_genomes.zip%22>`__

Then use kitsune kopt command below

\**Please be aware that this comand will use big computational resources
when large number of genomes and/or large genome size are used as the
input.

.. code:: bash

   kitsune kopt genome_list -kl 15 --canonical --fast -t 2 -o output.txt

.. |PyPI version| image:: https://badge.fury.io/py/kitsune.svg
   :target: https://badge.fury.io/py/kitsune
.. |Upload Python Package| image:: https://github.com/natapol/kitsune/workflows/Upload%20Python%20Package/badge.svg
