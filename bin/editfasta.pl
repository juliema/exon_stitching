#!/usr/bin/env perl

################################
#
#  This script edits the assembly atram output files. It removes wonky characters in the names of the genes
#        and adds a unique number to the end of each one. output is *.ed.fasta
#  This script also takes the file with all of the reference genes in it and makes a list of them and prints each
#        to its own file.  *.reference.fasta
#
################################


use strict;
use warnings;

my $pwd=shift;
my $query=shift;
#print "$pwd\n";

## get list of fasta  files in this directory
#`find -name "*best.fasta" | sed 's/^\\.\\///' >files`; # Use this version for very large numbers of files; sed patterns must be double-escaped
`find $pwd -empty -type f -delete`;
`find $pwd -name  "*.fasta"  >files`; # Use this version for very large numbers of files; sed patterns must be double-escaped

## open each file, remove funky characters and add a unique number to the end of each contig name. 
open FH, "<files"; 
while (<FH>) { # print;
	if (/(\S+).fasta/) {
		my $file=$1; 
		#print "$file\n";
		open FH2, "<$file.fasta";
		open OUT2, ">$file.ed.fasta";
		my $num=0;
		while (<FH2>) { #print;
			if (/^>(.*)$/) {	
				my $contig = $1;
				$num++;
 				$contig =~ s/[-,=@+\[\]:!]//g;
 				$contig =~ s/\s+//g;
				print OUT2 ">$contig\_$num\n";
			}
			elsif (/^\S+/ && ! /Command|Hostname|completed/) { print OUT2; }
		}
	}
}

open GENE_LIST, ">Gene_list.txt";
open QUERY, "<$query";
while (<QUERY>) {
	if (/^>(.*)/) { 	
		my $gene=$1; 
		$gene=~ s/\s+/_/g; 
		open GENE, ">$gene.reference.fasta";
		#print "$gene.fasta\n";
		print GENE_LIST "$gene.reference.fasta\n";
		print GENE ">$gene\n";
	}
	else { print GENE; }
}

#`rm files`;
