kitsune
--------

1. Kitsune
+
Kitsune is a toolkit for evaluate the lenght of optimal k-mer in dataset.

2. Basics
+
Kitsune provides many methods to evaluate the amount of information by given k-mer.
This quantitative measurement could be use to calculate the amount of need k-mer

3. Installation
+
--
a. install it throught pip
+
.................
$pip install kitsune
.................
__
4. Usage
+
--
a. Calculate k-mer from CRE, ACF, and OFC value
+
.................
$kitsune cre genome_fasta/* -ks 5 -ke 10
$kitsune acf genome_fasta/* -ks 5 -ke 10
$kitsune ocf genome_fasta/* -ks 5 -ke 10
.................

b. calculate deistance matrix
+
.................
$kitsune matrix genome_fasta/* -ks 5 -ke 10
.................
--