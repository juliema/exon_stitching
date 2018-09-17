# exon_stitching

This program will find and stitch together exons from targeted assemblies using amino acid targets and DNA assemblies.

### Usage: bash exon_stitching.sh  list_of_taxa.txt  AminoAcidFile.fasta  path_to_target_assemblies/ 


## Needs: 

1. a text file of all your taxon names
2. a fasta file of all the genes in amino acids 
3. the program exonerate -- https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide

The taxon and gene names must correspond to the assembly file names.  For example the assemblies from atram2.0 will produce a fasta file with the gene name followed by the taxon name. The taxon names and the gene names in the amino acid file must correspond. 

## Overview

This program will group all of the taxa for each gene together. It will use the program exonerate and your amino acid reference file to find exon positions in the assemblies. It will then stitch those exons together.  If there is a missing exon then the script will add groups of 3 NNNs in the missing places so that in the end, there will be 1 file for each gene with all of the exons for each taxon and with NNNS in the missing pieces. Therefore all of the exon lengths for each taxa are roughly the same, barring indels and ready for alignment steps.










