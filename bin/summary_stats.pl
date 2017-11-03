#!/usr/bin/env perl

use strict;
use warnings;

my $cwd = shift;

my @taxarray=();
my %taxhash=();
my %full=();
my %ninefive=();
my %ninezero=();
my %eightzero=();
my %sevenzero=();
my %fivezero=();
my %ten=();
my %lessten=();

open OUT, ">Summary_stats.per.taxon.csv";
print OUT "Taxon,NumGenes,FullExons,95%,90%,80%,70%,50%,10%\n";
print "Summary Stats Per Taxon\n";
print "Taxon,NumGenes,FullExons,95%,90%,80%,70%,50%,10%\n";
open FH, "<$cwd\/Summary_stats.per.gene.csv";
while (<FH>) {
	if (/\S+?,(\S+?),(\d+?),(\d+)$/) {
		my $tax=$1;
		#$tax =~ s/(\S+).ed.fasta/$1/;
		my $querylength=$2;
		my $targetlength=$3;
		my $proportion = $targetlength / $querylength;
		#print "$tax\t$proportion\n";
		#print "Proportion $proportion  Query length $querylength  Target Length  $targetlength  taxon $tax\n";
		if (! exists $taxhash{$tax}) {
			$taxhash{$tax}=1;
			$full{$tax}=0;
			$ninefive{$tax}=0;
			$ninezero{$tax}=0;
			$eightzero{$tax}=0;
			$sevenzero{$tax}=0;
			$fivezero{$tax}=0;
			$ten{$tax}=0;
			$lessten{$tax}=0;
			push @taxarray, $tax;
		}
		else {$taxhash{$tax}++;}
		if ($proportion == 1    ) { $full{$tax}++;     }
		elsif ($proportion >= 0.95 ) { $ninefive{$tax}++; }
                elsif ($proportion >= 0.90 ) { $ninezero{$tax}++; }
                elsif ($proportion >= 0.80 ) { $eightzero{$tax}++;}
                elsif ($proportion >= 0.70 ) { $sevenzero{$tax}++;}
                elsif ($proportion >= 0.50 ) { $fivezero{$tax}++; }
                elsif ($proportion >= 0.10 ) { $ten{$tax}++;      }   
                elsif ($proportion <  0.10 ) { $lessten{$tax}++;  }
	}
}
my $tot95=0;
my $tot90=0;
my $tot80=0;
my $tot70=0;
my $tot50=0;
my $tot10=0;
my $rest=0;

for my $tax (@taxarray) {
	#print "$taxhash{$tax} = number of genes\n";
	$tot95 = (  $ninefive{$tax}  + $full{$tax});
        $tot90 = (  $ninezero{$tax}  + $tot95);
        $tot80 = (  $eightzero{$tax} + $tot90);
        $tot70 = (  $sevenzero{$tax} + $tot80);
        $tot50 = (  $fivezero{$tax}  + $tot70);
        $tot10 = (  $ten{$tax}       + $tot50);
        $rest =  (  $lessten{$tax}   + $tot10);
	print OUT  "$tax,$taxhash{$tax},$full{$tax},$tot95,$tot90,$tot80,$tot70,$tot50,$tot10\n";
	print  "$tax,$taxhash{$tax},$full{$tax},$tot95,$tot90,$tot80,$tot70,$tot50,$tot10\n";
}	

