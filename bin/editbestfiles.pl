#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;
## get list of best files in this directory
#`find -name "*best.fasta" | sed 's/^\\.\\///' >files`; # Use this version for very large numbers of files; sed patterns must be double-escaped

`find $pwd -name  "*best.fasta"  >files`; # Use this version for very large numbers of files; sed patterns must be double-escaped

## open each file, remove funky characters and add a unique number to the end of ecach contig name. 
open FH, "<$pwd\/files";
while (<FH>) { # print;
	if (/(\S+).best.fasta/) {
		my $file=$1; print "$file\n";
		open FH2, "<$file.best.fasta";
		open OUT2, ">$file.best.ed.fasta";
		my $num=0;
		while (<FH2>) { #print;
			if (/^>(\S+)$/) {	
				my $contig = $1;
				$num++;
 				$contig =~ s/[-,=@+\[\]:!]//g;
				print OUT2 ">$contig\_$num\n";
			}
			elsif (/^\S+/ && ! /Command|Hostname|completed/) { print OUT2; }
		}
	}
}

