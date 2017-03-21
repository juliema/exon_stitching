#!/usr/bin/env perl

use strict;
use warnings;

my $overlap = 10;
my $file = 'results.stats.OVERLAP';
#PHUM615880.results.stats.OVERLAP.10.csv

`ls -l *.$file.$overlap.csv >files`;

#`ls -l *.results.stats.OVERLAP.$overlap.csv >files`;

#open OUT, ">SUMMARY.OVERLAP.csv";

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

open FH, "<files";
while (<FH>) {
	if (/(\S+).$file.$overlap.csv/) {
		my $gene=$1;
		$fullgenehash{$gene}=0;
		$partgenehash{$gene}=0;
		push @genearray, $gene;
		$totalgenes++;
		open FH1, "<$gene.$file.$overlap.csv";
		while (<FH1>) {
			if (/Allowing\s+overlap\s+(\S+)/) {
				 $overlap = $1;
			}
			if (/^(\S+?),/  && ! /Statistics/ && ! /Inputfile/ && ! /There\s+are/ && ! /Library/ && ! /Allowing/) {
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
				else { $libraryhash{$lib}++;
					#print "$gene\t$lib\n"; 
				}
			}
			#Aamic,3,206,FALSE,NA,trinity,3,1,1,trinity,205,0,205,
			#Hbarb,2,206,TRUE,trinity,,,,,,,,,2.c0_g2_i2_len=2115_path=1624039_17144053_17285455_460556144_4694145210_4760211764_5307765765_530876613
			if (/^(\S+?),\S+?,(\S+?),TRUE/) {
			#if (/^(Aamic),\S+?,(\S+?),TRUE/) {
				my $lib=$1;
				my $genelength=$2;
				#print "\t$lib FULL GENE\n";
				$libfullgenehash{$lib}++;
				$fullgenehash{$gene}++;  #print "full gene hash $gene $fullgenehash{$gene}\n";
			}						
			#Pnbad,3,206,FALSE,NA,trinity,3,1,1,trinity,156,0,156,3.c0_g1_i1_len=1095_path=107301094_35
			         #Aamic,  3 ,  206, FALSE,  NA,trinity,3,  1,   1,  trinity  ,205   ,0  ,205,
				#Library,NumContigs,GeneLength,(T/F),FG_Assembler,Assembler,LastContigNumber,TotalNumberofContigs,ContigstoKeep,Assembler,CombinedContigLength,TotalOverlap,Beginning,End,Beginning,End,ContigName
			elsif (/^(\S+?),\S+?,(\S+?),FALSE,\S+?,\S*?,\S*?,\S*?,(\S*?),(\S*?),\S*?,(\S*?),\S*?,/ && ! /Statistics/ && ! /Inputfile/ && ! /There\s+are/ && ! /Library/) {
				 # Hieur, 3,    191, FALSE,NA,trinity,3,   1,    1,   trinity, 0,   167,  24,191,CONTIGS,3.c1_g1_i1_len2976_path102975_1
			#elsif (/^(Aamic),\S+?,(\S+?),FALSE,\S+?,\S*?,\S*?,\S*?,(\S*?),(\S*?),(\S*?),\S*?,\S*?,/ && ! /Statistics/ && ! /Inputfile/ && ! /There\s+are/ && ! /Library/) {
				$partgenehash{$gene}++;
				my $lib=$1;
				#print "\t$lib partial gene\n";
				my $genelength=$2;
				my $numcontigs=$3;
				my $contassembler=$4;
				my $contiglength=$5;	
				my $proportion = $contiglength/$genelength;
				#print "$lib\t$genelength\t$contiglength\t$proportion\n";
				if ($proportion >= 0.95 ) { $ninefivehash{$lib}++; }
				if ($proportion >= 0.90 ) { $ninezerohash{$lib}++; } 
				if ($proportion >= 0.80 ) { $eightzerohash{$lib}++;} 
				if ($proportion >= 0.70 ) { $sevenzerohash{$lib}++;} 
				if ($proportion >= 0.50 ) { $fivezerohash{$lib}++; } 
				if ($proportion >= 0.10 ) { $onezerohash{$lib}++;  } 
				if ($proportion < 0.10 )  { $lesstenhash{$lib}++; }
				}		
		}
	}
}

print  "Number\tLibrary\tOverlap\tNumGenes\tNumFullGenes\tNum95\tNum90\tNum80\tNum70\tNum50\tNum10\tLessthan10\n";
my $count=0;
for my $lib (@libarray) {
	$count++;
	#print "$lib $libraryhash{$lib}\t$libfullgenehash{$lib}\t$ninefivehash{$lib}\t$ninezerohash{$lib}\n";
        my $tot95 = ($ninefivehash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot90 = ($ninezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot80 = ($eightzerohash{$lib}+ $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot70 = ($sevenzerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot50 = ($fivezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
        my $tot10 = ($onezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
	my $rest =  ($lesstenhash{$lib} + $onezerohash{$lib} + $libfullgenehash{$lib})/$libraryhash{$lib};
	print "$count\t$lib\t$overlap\t$libraryhash{$lib}\t$libfullgenehash{$lib}\t$tot95\t$tot90\t$tot80\t$tot70\t$tot50\t$tot10\t$rest\n";
}


