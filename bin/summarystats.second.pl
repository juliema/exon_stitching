#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;


my $overlap = 0;
my $file = 'exonerate2.stats.OVERLAP';

`ls -l $pwd\/*.$file.$overlap.csv >files`;
my $totalgenes=0;
my %libraryhash=();
my @genearray=();
my @libarray=();
my %partgenehash=();
my %libfullgenehash=();
my %ninefivehash=();
my %ninezerohash=();
my %eightzerohash=();
my %sevenzerohash=();
my %fivezerohash=();
my %onezerohash=();
my %lesstenhash=();
my %fullgenehash=();
my $countlib=0;
my %genecontighash=();

open FH, "<files";
while (<FH>) {  #print;
	if (/$pwd\/(\S+).$file.$overlap.csv/) {
		my $gene=$1;
		#print"$gene\n";
		$fullgenehash{$gene}=0;
		$partgenehash{$gene}=0;
		push @genearray, $gene;
		$totalgenes++;
		my $numlines=0;
		my $numfoundlines=0;
		open FH1, "<$pwd\/$gene.$file.$overlap.csv";
		while (<FH1>) { #print;
			if (/There\s+are\s+(\d+)\s+Libraries\s+in\s+(\S+)/) {
				my $numlib=$1;
				my $gene=$2;
				$genecontighash{$gene}=$numlib;
			}
			if (/Allowing\s+overlap\s+(\S+)/) {
				 $overlap = $1;
			}
			if (/^(\S+?),/  && ! /Statistics/ && ! /Inputfile/ && ! /There\s+are/ && ! /Library/ && ! /Allowing/) {
				$numlines++;
				my $lib = $1; # print "$lib\n";
				if (! exists $libraryhash{$lib}) {
					$libraryhash{$lib}=1;
					$countlib++;
					push @libarray, $lib;
					$ninefivehash{$lib}=0;
					$ninezerohash{$lib}=0;
					$eightzerohash{$lib}=0;
					$sevenzerohash{$lib}=0;
					$fivezerohash{$lib}=0;
					$onezerohash{$lib}=0;	
					$lesstenhash{$lib}=0;
					$libfullgenehash{$lib}=0;
					#print "$gene\t$lib\n";	
				}
				else { 
					$libraryhash{$lib}++;
					#print "$gene\t$lib\n"; 
				}
			}
			if (/^(\S+?),\S+?,(\S+?),TRUE/) {
				my $lib=$1;
				my $genelength=$2;
				$numfoundlines++;
				#print "\t$lib FULL GENE\n";
				$libfullgenehash{$lib}++;
				$fullgenehash{$gene}++;  #print "full gene hash $gene $fullgenehash{$gene}\n";
			}						
			my $checkflag = 0;
#			     PdhumHO,  3,    907,FALSE, NA,trinity,3,    0,   trinity,0,   0,  CONTIGS,
			if (/^(\S+?),\S+?,(\S+?),FALSE,\S+?,\S*?,\S*?,(\S*?),(\S*?),\S*?,(\S*?),\S*?,\S*?,/ && ! /Statistics/ && ! /Inputfile/ && ! /There\s+are/ && ! /Library/) {
				$partgenehash{$gene}++;
				$numfoundlines++;
				my $lib=$1;
				my $genelength=$2;
				my $numcontigs=$3;
				my $contassembler=$4;
				my $contiglength=$5;	
				my $proportion = $contiglength/$genelength;
				if ($proportion >= 0.95 ) { $ninefivehash{$lib}++;  $checkflag=1;}
				if ($proportion >= 0.90 ) { $ninezerohash{$lib}++;  $checkflag=1;} 
				if ($proportion >= 0.80 ) { $eightzerohash{$lib}++; $checkflag=1;} 
				if ($proportion >= 0.70 ) { $sevenzerohash{$lib}++; $checkflag=1;} 
				if ($proportion >= 0.50 ) { $fivezerohash{$lib}++;  $checkflag=1;} 
				if ($proportion >= 0.10 ) { $onezerohash{$lib}++;   $checkflag=1;} 
				if ($proportion < 0.10 )  { $lesstenhash{$lib}++;   $checkflag=1;}
				if ($checkflag == 0 ){ print "$lib\t$gene\t CHECKFLAG = 0, $genelength,  $contiglength, $proportion\n" };
			}		
		}
		my $diff = $numlines-$numfoundlines;
		if ($diff != 0) {
			print "PROBLEM\t$gene\t$numlines\t$numfoundlines\n";
		}
	}
}

print  "Number\tLibrary\tOverlap\tNumGenes\tNumFullGenes\tNum95\tNum90\tNum80\tNum70\tNum50\tNum10\tLessthan10\n";

my $count=0;
for my $lib (@libarray) {
	$count++;
	my $fullprop = $libfullgenehash{$lib}/$libraryhash{$lib};
        my $tot95 = ($ninefivehash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot90 = ($ninezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot80 = ($eightzerohash{$lib}+ $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot70 = ($sevenzerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot50 = ($fivezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot10 = ($onezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
	my $rest =  ($lesstenhash{$lib} + $onezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
	print "$count\t$lib\t$overlap\t$libraryhash{$lib}\t$fullprop\t$tot95\t$tot90\t$tot80\t$tot70\t$tot50\t$tot10\t$rest\n";
}


