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

### OUTPUTFILES >$gene.stitched_exons.fasta -- go directly to alignment
##### this prints out the stitched file and the summary file..

my $overlap = shift;
my $debug = 'debug';
open SUMMARY, ">Summary_stats.per.gene.csv";
print SUMMARY "Locus,Taxon,Query_Length,Target_Length\n";
##### go through and determine the number of loci then for each locus... do

## look for all the files with .contig_list in them..... 
`ls -l *.contig_list.csv >contig_files.txt`;
open FILES, "<contig_files.txt";
while (<FILES>) {
	if (/(\S+).overlap.*.contig_list.csv/) {
		my $gene=$1;  print "$gene\t$overlap\n $gene.overlap.$overlap.contig_list.csv\n\n\n\n\n\n\n";
		### print out a fasta file for each gene
		open FH, "<$gene.overlap.$overlap.contig_list.csv";	
		open OUT1, ">$gene.stitched_exons.fasta";
		while (<FH>) { #print;
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
					my $contig = $array[9];  if ($debug eq 'debug') { print "line 46 contig is $contig\n\n";}
					chomp $contig; $contig =~ s/,//g; $contig =~ s/\|//g;
					$contig = " " . $contig;   if ($debug eq 'debug') { print "line 48 contig is $contig\n\n";}
					my $seqflag=0;
					my $sequence=();
					open FASTA, "<$gene.exons.fasta";
					print OUT1 ">$infile.$contig.$numcontigs\n";  if ($debug eq 'debug') { print "printing $infile.$contig.$numcontigs  to outfile line 51\n"; }
					while (<FASTA>) {
						my $line=$_;
						#$line = " " . $line;
						### remove interleavedness
						chomp $line; #print "$line\n";
						if ($line =~ m/^>/) { #print "This should be the > line $line\n"; 
							$seqflag=0; 
							if ($debug eq 'debug') { print "line is $line  and infile is $infile, contig, $contig\n"; }
							if ($line =~ /$infile/) {
								if ($debug eq 'debug') {  print "line is $line  matched  $infile\n"; }
								my $linecontig = $line;
								$linecontig =~ s/>$infile\,(.*),\d+$/$1/g;  if ($debug eq 'debug') {  print "LINE 64 this should just be the contig $linecontig\n"; }
								$linecontig =~ s/,//g; $linecontig =~ s/\|//g; $linecontig = " " . $linecontig;  if ($debug eq 'debug') {  print "LINE 66 this should just be the contig $linecontig\n"; }
								if ($linecontig =~ /$contig/) {
									if ($debug eq 'debug') { print "THE infile is $infile  AND THE CONTIG IS $contig   MATCHED in this:\n$line\n\n";}
									$seqflag =1; if ($debug eq 'debug') { print "LINE 63 ONE CONTIG LOOP >$infile.$contig.$numcontigs\n\n";} 
								}
							}
						}
						elsif ($seqflag == 1 ) { $sequence = $sequence.$line;} #  print "adding line $line to sequence\n";}
					}
					close FASTA;
					if ($debug eq 'debug') { print "start $start end $end length $querylength\n";}
					###### PRINTS OUT NNNs AT THE BEGINNING OF THE GENE IF THE CONTIG DOES NOT START AT THE BEGINNING
					if ($start > 0) { for (0..$start) { print OUT1 "NNN"; } }
					print OUT1 "$sequence"; 
					##### PRINTS OUT NNNs AT THE END OF THE SEQUENCE IF CONTIG DOES NOT COVER THE FULL LENGTH
					if ($end < $querylength) { for ($end .. ($querylength-1)) { print OUT1 "NNN";  } }
					if ($debug eq 'debug') { my $num = ($querylength-1) - $end; print "printed $num groups of three NNNs\n";} 
					print OUT1 "\n";
				}
				######### ELSE IF MORE THAN ONE CONTIG STITCH THEM TOGETHER WITH NNNS IN THE MIDDLE 
				else {   if ($debug eq 'debug') { print "$infile more than one contig, stitching together numcontigs == $numcontigs\n";}
					my $count=0;
					my $contigstart=($numcontigs*2) + 7;
					my $gapstart=0;
					#start at 8 and go to 8+ number of contigs
					for ($contigstart..($contigstart + $numcontigs-1)) {
						my $contignumber=$_; if ($debug eq 'debug') { print "contig number is $contignumber\n";}
						my $start=$array[7+$count];
						my $end=$array[8+$count];			
						my $contig = $array[$contignumber];
						$contig = " " . $contig; $contig =~ s/\|//g;   if ($debug eq 'debug') { print "contig is $contig\n\n";}
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
						

									if ($line =~ m/^>/) { #print "This should be the > line $line\n"; 
										$seqflag=0; 
										if ($debug eq 'debug') { print "line is $line  and infile is $infile, contig, $contig\n"; }
										if ($line =~ /$infile/) {
											if ($debug eq 'debug') {  print "line is $line  matched  $infile\n"; }
											my $linecontig = $line;
											$linecontig =~ s/>$infile\,(.*),\d+$/$1/g;  if ($debug eq 'debug') {  print "this should just be the contig $linecontig\n"; }
											$linecontig = " " . $linecontig; $linecontig =~ s/\|//g; 
											if ($linecontig =~ /$contig/) {
												if ($debug eq 'debug') { print "THE infile is $infile  AND THE CONTIG IS $contig   MATCHED in this:\n$line\n\n";}
												$seqflag =1; if ($debug eq 'debug') { print "LINE 63 ONE CONTIG LOOP >$infile.$contig.$numcontigs\n\n";} 
												if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
												print OUT1 ">$infile.$contig.$numcontigs\n";  if ($debug eq 'debug') { print "printing name to outfile like 102\n"; }
											}
										}
									}
						




#									if ($line =~ /^>/) { 
#									$seqflag=0;
#									if ($line =~ /$infile/ && $line =~/$contig/) {
#										$seqflag =1;  #print OUT1;
#										##HERE 
#										print OUT1 ">$infile.$contig.$numcontigs\n";  if ($debug eq 'debug') { print "printing name to outfile like 102\n"; }
#										if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
#									}
#								}
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


								if ($line =~ m/^>/) { #print "This should be the > line $line\n"; 
									$seqflag=0; 
									if ($debug eq 'debug') { print "line is $line  and infile is $infile, contig, $contig\n"; }
									if ($line =~ /$infile/) {
										if ($debug eq 'debug') {  print "line is $line  matched  $infile\n"; }
										my $linecontig = $line;
										$linecontig =~ s/>$infile\,(.*),\d+$/$1/g;  if ($debug eq 'debug') {  print "this should just be the contig $linecontig\n"; }
										$linecontig = " " . $linecontig; $linecontig =~ s/\|//g; 
										if ($linecontig =~ /$contig/) {
											if ($debug eq 'debug') { print "THE infile is $infile  AND THE CONTIG IS $contig   MATCHED in this:\n$line\n\n";}
											$seqflag =1; if ($debug eq 'debug') { print "LINE 63 ONE CONTIG LOOP >$infile.$contig.$numcontigs\n\n";} 
										}		
									}
								}


#								if ($line =~ /^>/) { 
#									$seqflag=0;
#									if ($line =~ m/$infile/ && $line =~ /$contig/) {
#										$seqflag =1;  #print OUT1; 
#										#print OUT1 ">$infile.$contig.$numcontigs\n";
#										if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
#									}
#								}
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



								if ($line =~ m/^>/) { #print "This should be the > line $line\n"; 
									$seqflag=0; 
									if ($debug eq 'debug') { print "line is $line  and infile is $infile, contig, $contig\n"; }
									if ($line =~ /$infile/) {
										if ($debug eq 'debug') {  print "line is $line  matched  $infile\n"; }
										my $linecontig = $line;
										$linecontig =~ s/>$infile\,(.*),\d+$/$1/g;  if ($debug eq 'debug') {  print "this should just be the contig $linecontig\n"; }
										$linecontig = " " . $linecontig; $linecontig =~ s/\|//g; 
										if ($linecontig =~ /$contig/) {
											if ($debug eq 'debug') { print "THE infile is $infile  AND THE CONTIG IS $contig   MATCHED in this:\n$line\n\n";}
											$seqflag =1; if ($debug eq 'debug') { print "LINE 63 ONE CONTIG LOOP >$infile.$contig.$numcontigs\n\n";} 
										}	
									}
								}









								#if ($line =~ /^>/) { 
								#	$seqflag=0;
								#	if ($line =~ /$infile/ && $line =~ /$contig/) {
								#		$seqflag =1;  #print OUT1; 
								#		#print OUT1 ">$infile.$contig.$numcontigs\n";
								#		if ($debug eq 'debug') { print ">$infile.$contig.$numcontigs\n";}
								#	}
								#}
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
