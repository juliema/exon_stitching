#!/usr/bin/env perl

use strict;
use warnings;

my $pwd=shift;
my $overlap = 0;

### This script is for only one assembler when assembling genes. Other scripts are for if there are two and we want to take the best assembly.
### Contigs need to be in order they appear in the gene. so sort by the beginning of the contig.
### This script will look for all files in the folder called *.sorted, these should be exonerate results.
### It will create a file for each input file called  $input.stats.csv with all of the information about the contigs from exonerate
### rewriting script so there is only one loop for each of the multiple contigs.  so if Trinity > 1 go to same loop and if Trinity and Abyss > 1 go through the same Trinity loop.

#open STATS, ">OUTPUT.csv"; my $date=localtime(time);
`ls -l $pwd\/*exonerate2.sorted.ed.csv  >files`;
open FILES, "<files"; 
while (<FILES>) {
	if (/(\S+).exonerate2/) {
		my $inputfile=$1; 
		my $statoutput= "$inputfile.exonerate2.stats.OVERLAP.$overlap.csv";
		open STATS, ">$statoutput"; my $date=localtime(time);
		print STATS "Statistics from aTRAM assemblies exonerate2  $date\n";
		print STATS "Inputfile\t$inputfile.csv\n";
		# print "Statistics from aTRAM assemblies on $date\n";
		#print "Inputfile\t$inputfile.exonerate2.sorted.ed.csv\n";
		print STATS "Allowing overlap $overlap\n";
		#print  "Allowing overlap $overlap\n";
		my @taxarray=(); my $counttax=0; my %taxhash=();
		###############  Get list of taxa from file.
		open FH, "<$inputfile.exonerate2.sorted.ed.csv";
		#### debugging 11.17.16
		while (<FH>) { if (/^(\S+?),/ && ! /Library/) {  my $tax=$1; if (! exists $taxhash{$tax}) { $taxhash{$tax}=1; push @taxarray, $tax; $counttax++; } } }
		print STATS "There are $counttax Libraries in $inputfile\n";	
		print "There are $counttax Libraries in $inputfile\n";
		close FH;	
		print STATS "Library,Number of Contigs,GeneLength,Full Gene (T/F),FG_Assembler,Assembler,LastContigNumber,TotalNumberofContigs,ContigstoKeep,Assembler,CombinedContigLength,TotalOverlap,Beginning,End,Beginning,End,ContigName\n";
		############### Now for each taxon in the file loop through and find the best contig(s).
		for my $tax (@taxarray) {
			print STATS "$tax,"; print "\n\n$tax\n\n";
			my @contigarray=(); my $numcontigs=0; my @bigarray=(); my %contighash=();
			############  count the number of contigs for each taxon, and push the line to array of arrays, and contig into an array.
			open FH1, "<$inputfile.exonerate2.sorted.ed.csv";
			#print "opening $inputfile.exonerate2.csv.ed.sorted";
			my $pos=0;
			while (<FH1>) {  #print;
				if (/$tax?,(\S+?,\S+?,\S+?,\S+?,\S+?),(.*)$/ && ! /Library/) { print "found $tax adding 1 to $numcontigs\n";
					my $line = $1; my $contig = $2; push @contigarray, $contig; $numcontigs++; my @linearray = split(/,/, $line);
					for my $each (@linearray) { push @{$bigarray[$pos]}, $each;}
					$pos++;	
				}	
			}
			#print "The number of contigs is $numcontigs\n";
			close FH1;
			print STATS "$numcontigs,$bigarray[0][0],";  ## printing the number of contigs and the length of the full gene. 
			#####  Get statistics on each contig and find the best one per taxon.
			my $flag=0; #my $fullgene=0;
			my ($beg , $end) = (0) x 2;
 			my $assembler = $bigarray[0][4]; 
			for (0..($numcontigs-1)) {
				my $cont=$_;
				my $cdslength=$bigarray[$cont][0];
				my $alcontlength=$bigarray[$cont][1];
				#### if there is a contig with the whole length of the gene take it.
				if ($cdslength == $alcontlength) { 
					$contighash{$tax}=$contigarray[$cont];
					$flag=1;
					$beg=$bigarray[$cont][2];
					$end=$bigarray[$cont][3];
				}
			}
			####  if we have a full gene print out done.
			if ($flag == 1 ) { 
				print "fullgene\t$assembler and BEG $beg END $end\n"; print STATS "TRUE,$assembler,,,,,,,,0,$beg,$end,CONTIGS,$contighash{$tax}";  
			}
			#### If there is not one contig that is the whole length, go to last iteration, how many contigs from each assembler?  and examine the length of the contigs.
			if ($flag == 0) {
				# These variables are needed from every loop
				my ($sum, $numlastcontig,  $numtokeep, $lastcontig, $totoverlap) = (0) x 5;
				my (@begarray,  @endarray, @keepcontigarray, %keepcontighashBeg, %keepcontighashEnd);
				print STATS "FALSE,NA,";  ## not a full gene so FALSE in stats file and will give assembler later.
				#print "Line 83 Not full GENE number of contigs is $numcontigs\n\n";
				#for (0..($numcontigs-1)) {
				#	my $pos=$_;
					#if ($bigarray[$pos][4] eq 'trinity') {  
				#	my $contig=$contigarray[$pos]; 
				#	print "LINE 88 position $pos CONTIG $contig\n";
					#$contig =~ s/^(\d)\S+/$1/;
					#if ($contig == $lastcontig) { $numlastcontig++; }
					#elsif ($contig > $lastcontig) {  $numlastcontig=1; $lastcontig=$contig; }
				#}	
				print STATS "$assembler,$numcontigs,";
				#print  "$assembler,num contigs $numcontigs\n";   #Number of contigs last iteration $numlastcontig\n";
				#########################################################
				##### number of contigs = 1  
				#########################################################
				if ($numcontigs == 1) { #print "line 105 num last contigs =1\n";
					$numtokeep=1;  
					for (0..($numcontigs-1)) {
						my $pos=$_;
						#if ($bigarray[$pos][4] eq 'trinity') {
						#my $contig=$contigarray[$pos];
						#$contig =~ s/^(\d)\S+/$1/;
						#if ($contig == $lastcontig)  {
							$sum=$bigarray[$pos][1];
							my $contig=$contigarray[$pos]; push @keepcontigarray, $contig;
							my $beg=$bigarray[$pos][2]; push @begarray, $beg; $keepcontighashBeg{$contig}=$beg;
							my $end=$bigarray[$pos][3]; push @endarray, $end; $keepcontighashEnd{$contig}=$end;
							print STATS "1,$assembler,0,$sum,$beg,$end,CONTIGS,$contig";
							#print  "1,$assembler,0,$sum,$beg,$end,CONTIGS,$contig";
						#}
					}
				}
				my (@lengtharray); 
				#########################################################
				##### Number of Contigs  > 1  
				#########################################################
				if ($numcontigs > 1) {
					my @comparecontigarray=();
					#print "Number of contigs is > 1 total numcontigs = $numcontigs\n";       #and last contig = $numlastcontig\n";	
					my $total=$numcontigs-1;
					for (0..($total)) {
						my $pos=$_;
						#print "line 124 for 0 to $total position $pos\n";
						#if ($bigarray[$pos][4] eq 'abyss') {
						my $contig=$contigarray[$pos];
						#$contig =~ s/^(\d)\S+/$1/;
						#$tax last abyss  contig  $lastAcontig  number of contigs $numlastAcontig\n";
						#if ($contig == $lastcontig) {
							#print "last Abyss contig $contig == $lastAcontig\n";
						my $beg=$bigarray[$pos][2];
						my $end=$bigarray[$pos][3];
						push @begarray, $beg;
						push @endarray, $end;
						#my $keeptig = join("_", $contig,$beg,$end);
						#my $keeptig = $contigarray[$pos];
						push @comparecontigarray, $contig;
						push @lengtharray, $bigarray[$pos][1];
						#print "line 137 POSITION $pos contig $contig assembler $assembler length $bigarray[$pos][1]\t beginning  $beg\t end $end\n";
						#}    
					}
					my $pos=0;
					my $next=1;
					$totoverlap=0;
					my $truenumcontigs = $numcontigs-1;
					#### DEBUGGING 3.29.15
					#for (0..$numcontigs-1) {
					for (0..$truenumcontigs) {
						#print "line 148 POSITION $pos NEXT $next\n";
						#### DEBUGGING 3.29.15
						#print "$next <= ($truenumcontigs))\n";   #print "$next is <= $numcontigs-1\n";
						if ($next <= ($truenumcontigs)) {    #print "$next is <= $numlastcontig-1\n";
						#	if ($begarray[$next] > $endarray[$pos]) { print "line 153 does not overlap $begarray[$next] > $endarray[$pos]\n";	
							#### ADDING IN OVERLAP STEP: 			
				 			my $nextbeg = $begarray[$next]; my $endpos = $endarray[$pos];
							my $nextend = $endarray[$next];  #print "next beginning $nextbeg \t $nextend\n";
							if ($nextbeg > ($endpos-$overlap) && ($nextend > $endpos)) { 
								### CALCULATE OVERLAP  
								#print "does not overlap including $overlap buffer\n";
								if ($nextbeg < $endpos) {
									my $this = $endpos - $nextbeg;
									$totoverlap = $totoverlap + $this;
								}
								#print "line 153 does not overlap next beg $nextbeg > current end $endpos - minus overlap  $overlap and is not totally wihtin next end $nextend > $endpos\n";
								### Does not overlap, stitch together
								if ($sum == 0 ) { $sum = $lengtharray[$pos]+$lengtharray[$next]; }
								elsif ($sum > 0 ) { $sum = $sum +$lengtharray[$next]; } 
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
							elsif ($begarray[$next] <= ($endarray[$pos] - $overlap)) { 
								#print "overlaps $begarray[$next] <= $endarray[$pos] + $overlap\n"; print "sum = $sum\n";
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
										#$Askipsecond = 1;
										### keep the first 
										$next++;
										#print "line 205 $sum >= $secondsum keep original contigs\n";
										for my $contig (@keepcontigarray) {} #print "this contig $contig\n"; }
									}
									elsif ($sum < $secondsum) { #print "$sum < $secondsum keep second contig\n";
										$sum = $secondsum; 
										#print "Asum != 0  && Asum < secondsum abyss sum = $Asum\n";
										my $cont=$comparecontigarray[$next]; #print "THIS ONE $cont\n";
										undef (@keepcontigarray);
										undef (%keepcontighashBeg);
										undef (%keepcontighashEnd);
										$totoverlap=0;
										# DEBUGGING 3.27.15
										push @keepcontigarray, $cont; 
										#push @keepcontigarray, $contigarray[$next]; 
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
					#### DEBUGGING 4.6.14.  adding in tot overlap number
					#print STATS "num to keep $numtokeep,assembler $assembler,sum $sum,overlap $totoverlap,";
					print STATS "$numtokeep,$assembler,$totoverlap,$sum,";
					for my $contig  (@keepcontigarray) {
			    			print STATS "$keepcontighashBeg{$contig},$keepcontighashEnd{$contig},";
			    			#print  "sum $sum keeping  $keepcontighashBeg{$contig},$keepcontighashEnd{$contig},$contig,";
					}
					print STATS "CONTIGS,";
					for my $contig  (@keepcontigarray) {
			    			print STATS "$contig,";
					}
				}
			}
			##################################################################
			#          END!!!
			##################################################################
			print STATS "\n";
			undef(@contigarray);
			undef(@bigarray);
		}
	}
}
