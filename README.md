# LungCancerData
Early TCGA Data with RNAseq read depth curves for 4 genes. 

Data from Matt Wilkerson, Oct. 2011. 

Repository contains data (various formats), together with some simple analyses done in Matlab.

## Directories:
* OriginalData:  files from Matt
  * counts.csv:  Main file of read counts. Top row is case ID, 1st column is chromosome location, remaining columns are read counts
  * exonsMarron.csv:  Edited version, to fix errors in original file exons.csv  (by deleting the lines with the 1st column = 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 54.  Look at P10, Case Number 38)
* MatlabScripts:  used to generate files in all other directories (note this directory sturcture is needed for these to run properly)
* DerivedData:  Files generated by Matlab Scripts, e.g. convenient summaries in several formats
  * .mat, Matlab format, single file containing data for all 4 genes
  * .xls, Excel Spreadsheet version, one file for each gene
  * .cvs, Comma Separated Values, one file for each gene
* GraphicalOutputs:  Mostly useful visualizations of the data, and other analysis results


