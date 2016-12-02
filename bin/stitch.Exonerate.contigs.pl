#!/usr/bin/env perl                                                                                                                                                                                                                    
#####
#  THIS IS THE SECOND ROUND STITCHING SCRIPT. THIS SHOULD 
#    STITCH THE EXONERATE RESULTS TOGETHER WITH NNNS IN 
#      IN THE MIDDLE BEGINNING AND END IF IT IS NOT THE FULL CONTIG.
##########################################################

### OUTPUTFILES >$gene.Stiched.Final.Contigs.fasta -- go directly to alignment

### THIS SHOULD BE ZERO
my $overlap = 0;

system "ls -l *.exonerate2.stats.OVERLAP.$overlap.csv  >files";
open FH, "<files";
while (<FH>) {
	if (/(\S+).exonerate2.stats.OVERLAP.$overlap.csv/) {
		my $gene=$1;
		open FH1, "<$gene.exonerate2.stats.OVERLAP.$overlap.csv";
		open OUT1, ">$gene.Stiched.Final.Contigs.fasta";
		while (<FH1>) {  
			#print;	
			#### GET THE WHOLE CONTIG IF HAS THE WHOLE GENE
			##### EDITS APRIL 3 ONLY GET THE OPEN READING FRAME
			#Btmac,3,313,TRUE,trinity,,,,,,,,0,0,313,3.c
#			             $assembler,,,,,,,,0,$beg,$end,$contighash{$tax}
			    #PdhumCA,1,958,TRUE,trinity,,,,,,,,0,0,958,CONTIGS,PdhumCA_PHUM614880_whole_gene_0_958
			if (/^(\S+?),.*TRUE,\S+,,,,,,,,0,\d+,\d+,CONTIGS,(.*)$/) {
				my $lib=$1;
				my $contig=$2;
				print OUT1 ">$lib.$gene.1\n";
				#print "$lib\tTRUE and contig is $contig\n";
				my $seqflag=0;
				open FASTA, "<$gene.exonerate2.ed.fasta";
				while (<FASTA>) {
					$line=$_;
				#	print;
					chomp $line;
					if ($line =~ /^>/) {  $seqflag=0; }
					if ($seqflag == 1 ) { print OUT1 $line;}  # print $line;}
					if ($line =~ m/>$gene,\S+,$contig/g) {
						$seqflag =1;
						#print "MATCH $contig\n";
						#print OUT1 ">$lib\n";
			
		#print OUT1 ">TRUE$lib\_$contig\n";
						#print  ">$lib\_$contig\n";
					}
				}
				print OUT1 "\n";
			}
			#### DOES NOT HAVE THE WHOLE CONTIG
			#Aamic,5,958,FALSE,NA,trinity,5,4,trinity,0,723,2,285,376,504,606,760,782,940,CONTIGS,Aamic_PHUM614880_ONE_CONTIG_2_285,Aamic_PHUM614880_ONE_CONTIG_376_504,Aamic_PHUM614880_ONE_CONTIG_606_760,Aamic_PHUM614880_ONE_CONTIG_782_940, 
			elsif (/(\S+)/ && ! /Statistics/ && ! /Inputfile/ && ! /Allowing/ && ! /There\s+are/ && ! /Library/) {
				print "LINE 54 not full contig\n";
				my $line=$1;
				#print "$line\n";
				my @array = split (/,/, $line);
				#### FIND OUT IF ONLY NEED ONE CONTIG
				#Hbarb,2,313,FALSE,NA,trinity,1,2,1,trinity,310,0,3,313,1c
				my $numcontigs = $array[7];
				print "LINE 76 number of contigs is $numcontigs\n";
				#Hieur,1,958,FALSE,NA,trinity,1,1,trinity,0,268,2,270,CONTIGS,Hieur_PHUM614880_ONE_CONTIG_2_270
				if ($numcontigs == 1) {
					my $lib = $array[0];
					print "line 61  lib = $lib\n num contigs == 1\n";
					my $start=$array[11];
					my $end = $array[12];
					my $length=$array[2];
					print "\n$gene\t$lib ONE CONTIG $numcontigs START $start END $end LENGTH $length\n";
					my $contig = $array[14];
					chomp $contig;
					print "$contig\n";
					my $seqflag=0;
					my $sequence=();
					open FASTA, "<$gene.exonerate2.ed.fasta";
					while (<FASTA>) {
						$line=$_;
						#print;
						chomp $line;
						if ($line =~ /^>/) { $seqflag=0; }
						if ($seqflag == 1 ) { $sequence = $sequence.$line;}
						if ($line =~ m/^>$gene,\S+,$contig/) {
#						if ($line =~ m/^>$lib/) {
#						if ($line =~ m/$contig/) {
							$seqflag =1; 
						#	print "ONE CONTIG MATCH $lib  $contig\n";
							 #print OUT1; 
							print OUT1 ">$lib.$gene.1\n";
						#	print OUT1 ">$lib\n";
						#print;  #"\n>$lib\_$contig\n";
						}
					}
					###### PRINTS OUT NNNs AT THE BEGINNING OF THE GENE IF THE CONTIG DOES NOT START AT THE BEGINNING
					if ($start > 0) { for (0..$start) { print OUT1 "NNN"; } }
					print OUT1 "$sequence"; #print "$sequence";
					##### PRINTS OUT NNNs AT THE END OF THE SEQUENCE IF CONTIG DOES NOT COVER THE FULL LENGTH
					print "end $end length $length\n";
					if ($end < $length) { for ($end .. ($length-1)) { print OUT1 "NNN";  } } 
					print OUT1 "\n";
				}
				######### ELSE IF MORE THAN ONE CONTIG STITCH THEM TOGETHER WITH NNNS IN THE MIDDLE 
				else {
					my $lib = $array[0];
					print "\nlibrary $lib  more than one contig line 99\n";
					my $length = $array[2];
					$count=0;
					$contigstart=12;
					my $gapstart=0;
					print "LINE 123\n\n should be the whole line\n $line\n";
					$line =~ s/^\S+CONTIGS,(\S+)/$1/g;
					print "LINE 125 this should be the string of contigs\n$line\n";
					my @newcontigarray= split(/,/, $line);
					#for my $cont (@newcontigarray) { print "line 111 contig $cont\n"; }
					for (0..($numcontigs-1)) {
						my $contignumber=$_;
						print "line 113 Num contigs $numcontigs position $contignumber\n";
						#Aamic,5,958,FALSE,NA,trinity,5,4,trinity,0,723,2,285,376,504,606,760,782,940,CONTIGS
						my $start=$array[11+$count];
						my $end=$array[12+$count];			
						#my $pos=14 + (($numcontigs*2) + ($contignumber-1));
						#print "Position $pos\n";
						# REWRITING FOR NEW CONTIG ARRAY
						my $contig = $newcontigarray[$contignumber];
#						my $contig=$array[$pos];
#						my $contig=$array[$pos];
						print "$start,$end,$contig\n";
						my $last=0;
						my $gapend=0;
						#my $gapstart=0;
						my $sequence='';
						#print "CONTIG start $start end $end contig $contig\n";
						##### FIRST CONTIG
						if ( $contignumber == 0) { 
							print "$lib FIRST\t$start\t$end\n";
							print "contig $gene\t$contig\n";
							$gapstart=$end+1;
							### REMOVE
							#print OUT1 "\n";
							### REMOVE
							open FASTA, "<$gene.exonerate2.ed.fasta";
							while (<FASTA>) {
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line; }
								if ($line =~ m/>$gene,\S+,$contig/g) {
									$seqflag =1;  #print OUT1; 
									print OUT1 ">$lib.$gene.$numcontigs\n";
									#print OUT1 ">FIRST$lib\_$contig\_$start\_$end\n";
									#print;  #"\n>$lib\_$contig\n";
								}
							}
							#### AGAIN PRINT OUT NNNs IF NECESSARY AT THE BEGINNING
							if ($start != 0) { 
								for (0..$start) { print OUT1 "NNN";} 
							}
							print OUT1 "$sequence";
							$nextstart=$end;
							close FASTA;
						}
						$sequence='';
						### LAST CONTIG
						if ($contignumber == ($numcontigs-1))  { 
							#print "$lib LAST\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							#print "contig $gene\t$contig\n";
							#print OUT1 "GAPSTART\n";
							for ($gapstart ..($start-1)) { print OUT1 "NNN"; i} 
							print "LAST SEQ\n";
							#print OUT1 "  $lib.$contig\_$start\_$end\n";
							open FASTA, "<$gene.exonerate2.ed.fasta";
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$gene,\S+,$contig/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									print "\nLAST CONTIG >$lib  MATCH $contig\n";
								}
							}
							print "end $end  length $length\n";	
							print OUT1 "$sequence";
							#########  PRINT OUT NNNs AT THE END IF SEQUENCE DOES NOT COVER FULL LENGTH
							if ($end != $length) { print "end $end does not equal length $length NFILL\n";
								for ($end .. ($length-1)) { print OUT1 "NNN";  } 
							}
							print OUT1 "\n";
							close FASTA;
						}
						$sequence='';
						#### MIDDLE CONTIGS
						if ($contignumber > 0 && $contignumber < ($numcontigs-1)) {
							print "\nMIDDLE CONTIG\n";
							print "$lib MIDDLE\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							print "$lib\t$contig\n";
							#print OUT1 "GAPSTART\n";
							##### PRINT NNNs IN BETWEEN CONTIGS
							for ($gapstart ..($start-1)) { print OUT1 "NNN";} 
							#print OUT1 "LAST SEQ\n";
							#print "MATCH  $lib.$contig\_$start\_$end\n";
							open FASTA, "<$gene.exonerate2.ed.fasta";
							$gapstart = $end + 1;
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$gene,\S+,$contig/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									#print;  #"\n>$lib\_$contig\n";
								}
							}	
							print OUT1 "$sequence";
						}
						$count=$count+2;
					}
				}
			}
		}
	}
}

