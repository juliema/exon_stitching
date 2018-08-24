#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;
my $overlap = shift;
#my $debug = 'debug';
my $debug = 'off';

### Contigs need to be in order they appear in the gene. so sort by the beginning of the contig.
### This script will look for all files in the folder called *.sorted, these should be exonerate results.
### It will create a file for each input file called  $input.stats.csv with all of the information about the contigs from exonerate
### rewriting script so there is only one loop for each of the multiple contigs.  so if Trinity > 1 go to same loop and if Trinity and Abyss > 1 go through the same Trinity loop.

`ls -l $pwd\/*results.sorted.csv  >files`;
open FILES, "<files"; 
while (<FILES>) {
	if (/(\S+).results.sorted.csv/) {
		my $inputfile=$1;
		my $gene = $inputfile;
		$gene =~ s/\S+\///g;  if ($debug eq 'debug') {print "gene $gene\n";}
		###############  Create output file for contig names and where to stitch
		my $statoutput= "$inputfile.overlap.$overlap.contig_list.csv";
		open STATS, ">$statoutput"; my $date=localtime(time);
		#print STATS "Statistics from exon stitching  $date\n";
		print STATS "Inputfile\t$inputfile.csv   ";
		print STATS "Allowing overlap $overlap\n";
		my @taxarray=(); my $counttax=0; my %taxhash=();
		###############  Get list of taxa from file.
		open FH, "<$inputfile.results.sorted.csv";
		while (<FH>) { if (/^\S+?,(\S+?),/) {  my $tax=$1; if (! exists $taxhash{$tax}) { $taxhash{$tax}=1; push @taxarray, $tax; $counttax++; } } }
		print STATS "There are $counttax Libraries in $inputfile\n";	
		print STATS "Gene,Taxon,Number of Contigs,GeneLength,Contigs to Keep,Total Overlap,Combined Exon Length,Beginning,End,Beginning,End,ContigName\n";
		if ($debug eq 'debug') {print "There are $counttax Libraries in $inputfile\n";}
		close FH;	
		############### Now for each taxon in the file loop through and find the best contig(s).
		for my $tax (@taxarray) {
			print STATS "$gene,$tax,"; if ($debug eq 'debug') {print  "\n\n$tax\n";}
			my @contigarray=(); my $numcontigs=0; my @bigarray=(); my %contighash=();
			###############  Count the number of contigs for each taxon, and push the line to array of arrays, and contig into an array.
			open FH1, "<$inputfile.results.sorted.csv";
			my $pos=0;
			while (<FH1>) { #print;
				if (/\S+?,$tax,(\S+?,\S+?,\S+?,\S+?),(.*)$/) { if ($debug eq 'debug') {print "found $tax adding 1 to $numcontigs\n";}
					my $line = $1; my $contig = $2; push @contigarray, $contig; $numcontigs++; my @linearray = split(/,/, $line);
					for my $each (@linearray) { push @{$bigarray[$pos]}, $each;}
					$pos++;	
				}	
			}
			if ($debug eq 'debug') {print "The number of contigs is $numcontigs\n"; print " Number of contigs $numcontigs, length of full gene $bigarray[0][0]\n"; }
			close FH1;
			print STATS "$numcontigs,$bigarray[0][0],";  ## printing the number of contigs and the length of the full gene. 
			###############  Get statistics on each contig and determine if there is a full gene or not.
			my $flag=0;
			my ($beg , $end, $querylength,$queryoverlap) = (0) x 4;
			for (0..($numcontigs-1)) {
				my $cont=$_;
				my $targetlength=$bigarray[$cont][0];  if ($debug eq 'debug') {print "target exon length $targetlength\n"; }
				$querylength=$bigarray[$cont][1];  if ($debug eq 'debug') {print "qurey exon length $querylength contig number $cont\n";}
				#### if there is a contig with the whole length of the gene take it.
				if ($targetlength == $querylength) { if ($debug eq 'debug') {print "Target length == the query length, full gene found\n";}
					$contighash{$tax}=$contigarray[$cont];
					$flag=1;
					$beg=$bigarray[$cont][2]; if ($debug eq 'debug') {print "beginning location $beg\n";}
					$end=$bigarray[$cont][3]; if ($debug eq 'debug') {print "ending location $end\n";}
				}
			}
			####  if we have a full gene print out done.
			if ($flag == 1 ) { 
				print STATS "1,$queryoverlap,$querylength,$beg,$end,$contighash{$tax}";  
			}
			#### If there is not one contig that is the whole length.
			if ($flag == 0) {
				my ($sum, $numlastcontig,  $numtokeep, $lastcontig, $totoverlap) = (0) x 5;
				my (@begarray,  @endarray, @keepcontigarray, %keepcontighashBeg, %keepcontighashEnd);
				my $total=$numcontigs-1; 
				#########################################################
				##### number of contigs = 1  
				#########################################################
				if ($numcontigs == 1) { if ($debug eq 'debug') {print " One contig not full gene\n";}
					$numtokeep=1;  
					for (0..($total)) {
						my $pos=$_;
						$sum=$bigarray[$pos][1];
						my $contig=$contigarray[$pos]; push @keepcontigarray, $contig;
						my $beg=$bigarray[$pos][2]; push @begarray, $beg; $keepcontighashBeg{$contig}=$beg;
						my $end=$bigarray[$pos][3]; push @endarray, $end; $keepcontighashEnd{$contig}=$end;
						print STATS "$numcontigs,$queryoverlap,$sum,$beg,$end,$contig";  if ($debug eq 'debug') {print "$tax one contig total $sum  beginning $beg end $end contig $contig\n";}
					}
				}
				my (@lengtharray); 
				#########################################################
				##### Number of Contigs  > 1  
				#########################################################
				if ($numcontigs > 1) {   if ($debug eq 'debug') {print "$tax number of contigs is >1 num contigs = $numcontigs\n";}
					my @comparecontigarray=();
					for (0..($total)) {
						my $pos=$_;
						my $contig=$contigarray[$pos]; if ($debug eq 'debug') {print "contig is $contig\t";}
						my $beg=$bigarray[$pos][2];    if ($debug eq 'debug') {print "beginning  is $beg\t";}
						my $end=$bigarray[$pos][3];    if ($debug eq 'debug') {print "end  is $end\n";}
						push @begarray, $beg; push @endarray, $end; push @comparecontigarray, $contig;
						push @lengtharray, $bigarray[$pos][1];  if ($debug eq 'debug') {print "length is $bigarray[$pos][1]\n";}
					}
					my $pos=0;
					my $next=1;
					for (0..$total) {
						if ($next <= ($total)) {   
				 			my $nextbeg = $begarray[$next]; if ($debug eq 'debug') {print "beginning  is $nextbeg\t";}
							my $endpos = $endarray[$pos];
							my $nextend = $endarray[$next]; if ($debug eq 'debug') {print "end  is $nextend\n";}  
							if ($nextbeg > ($endpos-$overlap) && ($nextend > $endpos)) { if ($debug eq 'debug') {print "beginign of next is greater than the end of the last minus overlap. and the next end is greater than the last end. Does not overlap.\n";}
									### calculate overlap -- might be able to simpifly this!	
								if ($nextbeg < $endpos) {
									my $this = $endpos - $nextbeg;
									$totoverlap = $totoverlap + $this; if ($debug eq 'debug') {print "overlap calculated is $totoverlap\n";}
								
								}
								### Does not overlap, stitch together
								if ($sum == 0 ) { $sum = $lengtharray[$pos]+$lengtharray[$next]; }
								elsif ($sum > 0 ) { $sum = $sum +$lengtharray[$next]; } 
								if ($debug eq 'debug') {print "sum of query is $sum\n";}
								#print "stitch together sum is $sum\n";
								my $contF=$comparecontigarray[$pos];
								my $contS=$comparecontigarray[$next];
								if (! exists $keepcontighashBeg{$contF} ) {
									push @keepcontigarray, $contF;
									$keepcontighashBeg{$contF}=$begarray[$pos];
									$keepcontighashEnd{$contF}=$endarray[$pos];
									$numtokeep=$numtokeep+1;			
								}
								if (! exists $keepcontighashBeg{$contS} ) {
									push @keepcontigarray, $contS;
									$keepcontighashBeg{$contS}=$begarray[$next]; 
									$keepcontighashEnd{$contS}=$endarray[$next];
									$numtokeep=$numtokeep+1;			
								}
								$pos=$next;
								$next++;
							 }
							### overlaps, find the longest;
							elsif ($nextbeg <= ($endpos - $overlap)) { 
								if ($sum == 0) { #print "There are no contigs kept yet $sum = 0\n";	
									my $firstsum=$lengtharray[$pos];
									my $secondsum=$lengtharray[$next];
									if ($firstsum >= $secondsum) { #print "$firstsum >= $secondsum keep the fist contig OVERLAP IS $totoverlap\n";
										$sum = $firstsum; 
										my $cont=$comparecontigarray[$pos];
										push @keepcontigarray,  $cont;   #$Acontigarray[$pos];
										$keepcontighashBeg{$cont} = $begarray[$pos];  #print "line 182  KEEP FIRST ADDING BEGINNING $begarray[$pos]\t";
										$keepcontighashEnd{$cont} = $endarray[$pos];  #print " line 183 KEEP FIRST ADDING END  $endarray[$pos]\n";
										$numtokeep=$numtokeep+1;
									}
									elsif ($firstsum < $secondsum) { #print "$firstsum < $secondsum  keep second contig\n";
										$sum = $secondsum;
										######## KEEP SECOND CONTIG	
										my $cont=$comparecontigarray[$next];
										push @keepcontigarray, $cont;   #$Acontigarray[$next];
										$keepcontighashBeg{$cont} = $begarray[$next];  #print " line 191 KEEP SECOND ADDING BEGINNING $begarray[$next]\t";
										$keepcontighashEnd{$cont} = $endarray[$next];  #print " line 192 KEEP SECOND  ADDING END $endarray[$next]\n";
										$numtokeep=$numtokeep+1;
										$pos=$next;
									}
									$next++;
								}
								elsif ($sum != 0) { #print "line 201 Already have contigs in memory\n";
									my $secondsum=$lengtharray[$next];
									if ($sum >= $secondsum) {
										$next++;
										for my $contig (@keepcontigarray) {} #print "this contig $contig\n"; }
									}
									elsif ($sum < $secondsum) { 
										$sum = $secondsum; 
										my $cont=$comparecontigarray[$next]; #print "THIS ONE $cont\n";
										undef (@keepcontigarray);
										undef (%keepcontighashBeg);
										undef (%keepcontighashEnd);
										$totoverlap=0;
										push @keepcontigarray, $cont; 
										$keepcontighashEnd{$cont} = $endarray[$next];  #print " line 218 KEEP SECOND  ADDING END $endarray[$next]\n";
										$keepcontighashBeg{$cont} = $begarray[$next];  #print  " line 217 KEEP SECOND ADDING BEGINNING $begarray[$next]\t";
										$numtokeep=1;
										$pos=$next;
										$next++;
									}
								}			    		
							}
						}
					}
					$numtokeep=scalar(@keepcontigarray); #print "line 232 $numtokeep contigs to keep\n"; print "sum = $sum\n";
					for (0..$numtokeep-1) { $pos=$_;} #print "line 227 KEEPING $pos $keepcontigarray[$pos]\n"; } 
					##################################################################
					#          FINAL PRINT
					##################################################################
					print STATS "$numtokeep,$totoverlap,$sum,";  if ($debug eq 'debug') {print "num to keep $numtokeep,total overlap $totoverlap, sum $sum,";}
					for my $contig  (@keepcontigarray) {
			    			print STATS "$keepcontighashBeg{$contig},$keepcontighashEnd{$contig},";  if ($debug eq 'debug') {print "$keepcontighashBeg{$contig},$keepcontighashEnd{$contig},";}
					}
					for my $contig  (@keepcontigarray) {
			    			print STATS "$contig,"; if ($debug eq 'debug') {print "$contig,";}
					}
				}
			}
			##################################################################
			#          END!!!
			##################################################################
			print STATS "\n"; if ($debug eq 'debug') {print "\n";}
			undef(@contigarray);
			undef(@bigarray);
		}
	}
}
