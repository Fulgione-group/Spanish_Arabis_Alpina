###vei
# Java programs 
###

Java programs to produce a SNP matrix file from a vcf, and to compute the density of pairwise differences across all sample pairs from the SNP matrix.
Note that the matrix has to include segregating and non-segregating sites.

These java programs assume that you have the /java/ folder from GitHub in your home directory, and that the directory structure is the same as in /java/. 
If not so, you will need to change some paths (e.g. to the libraries in ./java/lib/)


To run these java programs, enter the ./java/projects directory and you will find:

###  matrix_alpina.command
Reads a VCF file and creates a SNP matrix filtered for coverage and quality. The SNP matrix is a lot smaller than a vcf file, which speeds up analyses.

### pi_matrix_alpina.command
Computes the density of pairwise differences between a focus sample and all others, based on the SNP matrix created by matrix_alpina.command
You should imagine the output of this as a single line in a matrix of pairwise differences across all sample pairs (n x n).
