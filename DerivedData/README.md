# Derived Data Files
Early TCGA Data with RNAseq read depth curves for 4 genes. 

Data from Matt Wilkerson, Oct. 2011. 

The 4 genes are:
* PIK3CA
* CDKN2A
* P10
* STK11

Files in this directory are intended for easy downloading and use

##  File formats
* All are matrices
* All column headers are Case Names  (TCGA IDs)
  * Note these are the same for all 4 genes
* All row headers are Genomic Coordinates of Base Pairs
* All main data entries are read counts
* Matlab .mat file has all three genes in the same file
  * Genes can be separated, and worked with, using read statements as in LungCancerViz.m in the MatlabScripts directory
* Remaining file types are separate for each gene
  *  Gene names are appended at the end of each filename 



