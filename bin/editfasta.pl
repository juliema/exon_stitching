#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;
my $query=shift;
#print "$pwd\n";

## get list of fasta  files in this directory
#`find -name "*best.fasta" | sed 's/^\\.\\///' >files`; # Use this version for very large numbers of files; sed patterns must be double-escaped

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
		open GENE, ">$gene.fasta";
		print GENE_LIST "$gene.fasta\n";
		print GENE ">$gene\n";
	}
	else { print GENE; }
}

