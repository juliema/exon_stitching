#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;

`ls -l $pwd/*.exonerate2.sorted.csv  >files`;

open FH, "<files";
while (<FH>) {
	if (/(\S+).exonerate2.sorted.csv/) {
		my $gene=$1;
		open FH1, "<$gene.exonerate2.sorted.csv";
		open OUT, ">$gene.exonerate2.sorted.ed.csv";
		while (<FH1>) {
			if (/(\d+?,\d+?),(\d+?,\d+?),(\S+?),((\S+_\S+?)_\S+)$/) {
				my $first=$1;
				my $info=$2;
				my $next=$3;
				my $contig=$4;
				my $tax=$5;
				print OUT "$tax,$first,$info,$next,$contig";
				$info =~s /,/_/g;
				print OUT "_$info\n";
			}			
		}
		open FH2, "<$gene.exonerate2.fasta";
		open OUT1, ">$gene.exonerate2.ed.fasta";
		while (<FH2>) {
			if (/^>(\S+),(\S+,\d+,\d+)$/) {
				my $gene=$1;
				my $contig=$2;
				$contig =~ s/,/_/g;
				print OUT1 ">$gene,$contig\n";
			}
			else {  print OUT1; }
		}
	}
}
