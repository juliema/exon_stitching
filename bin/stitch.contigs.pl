#!/usr/bin/env perl   

use strict;
use warnings;

##########################################################
#
#  	THIS SCRIPT  STITCHES CONTIGS TOGETHER. THIS SHOULD 
#   	STITCH THE EXONERATE RESULTS TOGETHER WITH NNNS IN 
#      	IN THE MIDDLE BEGINNING AND END IF IT IS NOT THE FULL CONTIG.
#
##########################################################

### OUTPUTFILES >$gene.Stiched.Final.Contigs.fasta -- go directly to alignment
##### this prints out the stitched file and the summary file..

my $overlap = shift;
my $debug = 'debug';
open SUMMARY, ">Summary_stats.csv";
print SUMMARY "Locus,Taxon,Query_Length,Target_Length\n";
##### go through and determine the number of loci then for each locus... do

## look for all the files with .contig_list in them..... 
`ls -l *.contig_list.csv >contig_files.txt`;
open FILES, "<contig_files.txt";
while (<FILES>) {
	if (/(\S+).overlap.*.contig_list.csv/) {
		my $gene=$1;
		### print out a fasta file for each gene
		open FH, "<$gene.overlap.$overlap.contig_list.csv";	
		open OUT1, ">$gene.stitched_exons.fasta";
		while (<FH>) { print;
			### this line look for the gene name at the beginning....
			if (/(\S+)/ && ! /Inputfile/ && ! /There\s+are/ && ! /Number\s+of\s+Contigs/) {
				my $line=$1;  if ($debug eq 'debug') { print "line from stats file is $line\n\n\n";}
				my @array = split (/,/, $line);
				#### get number of contigs
				my $numcontigs = $array[4];   if ($debug eq 'debug') { print "LINE 38 num contigs is $numcontigs\n";}
				my $querylength=$array[3];
				my $combinedexonlength=$array[6];
				my $infile = $array[1]; if ($debug eq 'debug') { print "LINE 43 infile from stats file is $infile\n\n";}
				print SUMMARY "$gene,$infile,$querylength,$combinedexonlength\n";
				if ($numcontigs == 1) {   if ($debug eq 'debug') { print "LINE 46 infile $infile one contig $numcontigs\n\n";} 
					my $start=$array[7];
					my $end = $array[8];
					my $contig = $array[9];  if ($debug eq 'debug') { print "contig is $contig\n\n";}
					chomp $contig; $contig =~ s/,//g;
					      my $seqflag=0;
					my $sequence=();
					open FASTA, "<$gene.exons.fasta";
					while (<FASTA>) {
						my $line=$_;
						#print "$line";
						### remove interleavedness
						chomp $line;
						if ($line =~ m/^>/) { 
							$seqflag=0; 
							if ($line =~ m/$infile,$contig/g) {
								if ($debug eq 'debug') { print "THE infile is $infile  AND THE CONTIG IS $contig   MATCHED in this:\n$line\n\n";}
								$seqflag =1; if ($debug eq 'debug') { print "LINE 63 ONE CONTIG LOOP >$infile.$contig.$numcontigs\n\n";} 
								print OUT1 ">$infile.$contig.$numcontigs\n";
							}
						}
						elsif ($seqflag == 1 ) { $sequence = $sequence.$line;} #  print "adding line $line to sequence\n";}
					}
					close FASTA;
					#if ($debug eq 'debug') { print "start $start end $end length $querylength\n";}
					###### PRINTS OUT NNNs AT THE BEGINNING OF THE GENE IF THE CONTIG DOES NOT START AT THE BEGINNING
					if ($start > 0) { for (0..$start) { print OUT1 "NNN"; } }
					print OUT1 "$sequence"; #print "SEQUENCE IS $sequence\n";
					##### PRINTS OUT NNNs AT THE END OF THE SEQUENCE IF CONTIG DOES NOT COVER THE FULL LENGTH
					if ($end < $querylength) { for ($end .. ($querylength-1)) { print OUT1 "NNN";  } }
					if ($debug eq 'debug') { my $num = ($querylength-1) - $end; print "printed $num groups of three NNNs\n";} 
					print OUT1 "\n";
				}
				######### ELSE IF MORE THAN ONE CONTIG STITCH THEM TOGETHER WITH NNNS IN THE MIDDLE 
				#$count=$count+2;
				else {   if ($debug eq 'debug') { print "$infile more than one contig, stitching together numcontigs == $numcontigs\n";}
					#print "\nlibrary $lib  more than one contig line 99\n";
					my $count=0;
					my $contigstart=($numcontigs*2) + 7;
					my $gapstart=0;
					#start at 8 and go to 8+ number of contigs
					for ($contigstart..($contigstart + $numcontigs-1)) {
						my $contignumber=$_; if ($debug eq 'debug') { print "contig number is $contignumber\n";}
						my $start=$array[7+$count];
						my $end=$array[8+$count];			
						my $contig = $array[$contignumber];
						my $last=0;
						my $gapend=0;
						my $sequence='';
						my $seqflag='';
						if ($debug eq 'debug') { print "contig number $contignumber, start = $start end = $end\n";}	
						### FIRST CONTIG
						if ( $contignumber == $contigstart) { if ($debug eq 'debug') { print "contig number equals contig start\n";}
							$gapstart=$end+1;
							open FASTA, "<$gene.exons.fasta";
							while (<FASTA>) {
								my $line=$_;
								chomp $line;
								if ($line =~ /^>/) { 
									$seqflag=0;
									if ($line =~ m/$infile,$contig/g) {
										$seqflag =1;  #print OUT1; 
										print OUT1 ">$infile.$contig.$numcontigs\n";
										if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
									}
								}
								elsif ($seqflag == 1 ) { $sequence = $sequence.$line; }
							}
							close FASTA;
							#### AGAIN PRINT OUT NNNs IF NECESSARY AT THE BEGINNING
							if ($start != 0) { 
								for (0..$start) { print OUT1 "NNN";} 
							}
							print OUT1 "$sequence";
							my $nextstart=$end;
							$gapstart =$end+1; 
						}
						### LAST CONTIG
						elsif ($contignumber == (($numcontigs-1) + $contigstart))  { if ($debug eq 'debug') { print "last contig $contignumber\n";}
							if ($debug eq 'debug') { print "start of gap is $gapstart and it goes until $start-1\n"; }	
							## make sure overlap is accounted for....
							for ($gapstart ..($start-1)) { print OUT1 "NNN";} 
							open FASTA, "<$gene.exons.fasta";
							while (<FASTA>) {
								my $line=$_;
								chomp $line;
								if ($line =~ /^>/) { 
									$seqflag=0;
									if ($line =~ m/$infile,$contig/g) {
										$seqflag =1;  #print OUT1; 
										#print OUT1 ">$infile.$contig.$numcontigs\n";
										if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
									}
								}
								elsif ($seqflag == 1 ) { $sequence = $sequence.$line; }
							}
							print OUT1 "$sequence";
							#########  PRINT OUT NNNs AT THE END IF SEQUENCE DOES NOT COVER FULL LENGTH
							if ($end < $querylength) { #print "end $end does not equal length $length NFILL\n";
								for ($end .. ($querylength-1)) { print OUT1 "NNN";  } 
							}
							print OUT1 "\n";
							close FASTA;
						}
						#### MIDDLE CONTIGS
						elsif ($contignumber > $contigstart && $contignumber < ($numcontigs-1 + $contigstart)) {  if ($debug eq 'debug') { print "middle contig $contignumber\n";}
							### here take into account overlap.
							if ($debug eq 'debug') { print "start of gap is $gapstart and it goes until $start-1\n"; }
							for ($gapstart ..($start-1)) { print OUT1 "NNN";} 
							open FASTA, "<$gene.exons.fasta";
							$gapstart = $end + 1;
							while (<FASTA>) {
								my $line=$_;
								chomp $line;
								if ($line =~ /^>/) { 
									$seqflag=0;
									if ($line =~ m/$infile,$contig/g) {
										$seqflag =1;  #print OUT1; 
										#print OUT1 ">$infile.$contig.$numcontigs\n";
										if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
									}
								}
								elsif ($seqflag == 1 ) { $sequence = $sequence.$line; }
							}
							print OUT1 "$sequence";
						}
						$count = $count + 2;
					}
				}
			}
		}
	}
}
