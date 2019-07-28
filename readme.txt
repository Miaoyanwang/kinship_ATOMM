#########################################################################################
Kinship_ATOMM is a C program that calculates genetic relatedness matrix among haploid individuals using both mutation and deletion polymorphisms.

The program takes the following as input:

1. genotypefile containing the genotype information for the studied individuals.
2. sizefile listing the number of individuals, the number of genetic variants, and the minor allele frequency threshold for the analysis.  

The problem produces two output files, here called ``kinshipfile'' and ``freqfile".
1. kinshipfile giving the genetic relatedness matrix for the studied individuals. 
2. freqfile listing the allele frequencies for the studied variants. 

Software accompaniment to:

"Two-way Mixed-Effects Methods for Joint Association Analyses Using Both Host and Pathogen Genomes. 
M. Wang, F. Roux, C. Bartoli, C. H.-Chauveau, C. Meyer, H. Lee, D. Roby, M. S. McPeek, and J. Bergelson Proc. Natl. Acad. Sci. Vol. 115 (24), E5440-E5449, 2018."

#########################################################################################
Installation:

(1) Download the kinship_ATOMM package. This package contains an example folder and the following source code:

kinship_ATOMM.cnrutil.cnrutil.h

The package also contains a pre-compiled binary executable file, called kinship_ATOMM.o, for MAC.

(2) If the pre-compiled binary file is not suitable for your use, please compile the program on your own machine. 
A C compiler should be used to compile the program. When using the gcc compiler, type the following in the terminal: 

gcc kinship_ATOMM.c -o kinship_ATOMM.o

(3) Run kinship_ATOMM

The program takes the following as input:

1. genotypefile containing the genotype information for the studied individuals.
2. sizefile listing the number of individuals, the number of genetic variants, and the minor allele frequency threshold for the analysis.  

The problem produces two output files, here called ``kinshipfile'' and ``freqfile".
1. kinshipfile giving the genetic relatedness matrix for the studied individuals. 
2. freqfile listing the allele frequencies for the studied variants.   

The format of the input and output will be described in the next section. To run the executable problem, type the following command:

./kinship_ATOMM.o -g genotypefile.txt -s sizefile.txt -k kinshipfile.txt -f freqfile.txt

The above command takes genotypefile.txt, sizefile.txt as inputs and generates kinshipfile.txt, freqfile.txt as outputs. 

###############################################################
Input:

1) Format for genotype file

This input file should consist of genotype information for all studied individuals. The file consists of (m+1) rows and (n+1) columns, where m is the number of genetic variants and n is the number of individuals. The first row of the file is a header that lists the individual ID. Starting form the second row, each row corresponds to a genetic variant. The first column of the file gives the genetic variant ID. Starting from the second column, each column corresponds to an individual. 

The genotype should be encoded in one of the followings:

A, T, C, G, -, and, N

Here "-" represents deletion, "N" represents missing data, and "A","T","C","G" represent the allele. 


Example: 
index BR01 BR02 BR03 BR04 BR06 BR18 BR19 BR20 BR21 BR14 BR16 BR13 BR15 BR07 BR23 BR24 BR25 BR08 BR09 BR10 BR11 BR12
1 - - - - - - - - - - - A - - - - G - - - - -
2 - - - - - - - - - - - G - - - - C - - - - -
(1)(2)...(n+1)

(1) genetic variant ID
(2) genotype for individual 1
(3) genotype for individual 2
...
(n+1) genotype for individual n

The genetic variant ID does not have to be consecutive. Please note that a blank space should be included between any two columns. 


2) Format for size file

This input should contain three numerical values in a row. The three numbers are, respectively, the number of individuals, the number of genetic variants, and the desired minor allele frequency (MAF) threshold. The number of individuals and the number of genetic variants should match with the count in the genotype file. 

Example:
22 144700 0.05
(1) (2) (3)

(1) the total number of individuals
(2) the total number of genetic variants
(3) the desired MAF threshold. Genetic variants with MAF below than this threshold will be removed from the analysis. 

###############################################################
Output:

1) Format for kinship file
This output contains the genetic relatedness matrix (GRM) for the studied individuals. The GRM is of size n-by-n, where n is the number of individuals in the input file. Each entry in the matrix corresponds to the genetic relatedness among the pair of individuals. 

Example:
 1.164182 -0.577799 ...-0.611849 -0.583317
-0.577799  1.111720 ... 0.649917  0.737269
-0.613349  0.601231 ... 0.527675  0.598031
-0.613604  0.602301 ... 0.526585  0.597657
	

2) Format for frequency file
This output contains the allele frequency for the studied variants. Each row gives the allele frequency at a genetic variant. 

Example:
index	 allel1_freq	 allel2_freq	 deletion_freq
1	0.045	0.045	0.909
2	0.045	0.045	0.909
(1)	(2)	(3)	(4)

(1) genetic variant ID
(2) frequency for the 1st allele 
(3) frequency for the 2nd allele 
(4) frequency for the deletion

###############################################################
Example:
The example folder contains example inputs and outputs. Input files are "genotype.txt", "size.txt". Output files are "freq.txt" and "kinship.txt"

To run the program using the example input files, type:

./kinship_ATOMM.o -g example/input/genotype.txt -s example/input/size.txt -f freq.txt -k kinship.txt

Compare the results to the example output files in the ``example/output'' folder. 