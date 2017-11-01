#!/bin/bash

########################################################################################
#Usage: sh exon_stitching.sh AminoAcidFile.fasta path_to_target_assemblies/ output_name
########################################################################################

# help menu
if [ "$1" = "-h" ] ; then
  echo "Usage: `basename $0` AminoAcidFile.fasta  path_to_target_assemblies/  overlap";
  exit 0
fi

echo amino acid file $1;
echo path to assemblies $2; 
myref=$1;
fastafilespath=$2;
myoverlap=$3;  
cwd=$(pwd);
mygene=${myref%.fasta} 
echo $mygene;
summaryfile='SummaryFile.txt';  
numgenes=0;


#step 1
#echo "step 1: edit fasta file names";
#echo targeting the following assemlies;
perl $cwd\/bin/editfasta.pl $fastafilespath;
#echo "done with step 1";

##step 2
#### check to see if it is just one fasta file or a list of genes.
if [[ $myref ==  *".fasta"* ]]; then
	numgenes=1;
	#echo "number of genes is $numgenes";
fi


##########################################################
#### List of genes now then loop through each gene for exonerate
########################################################

echo "step 2: running exonerate";
if [[ $numgenes == 1 ]]; then 

	echo "number of genes ==  1";
	for taxon in $fastafilespath*.ed.fasta; do
		#echo $lib $fastafilespath;
		#tax=${lib#"$fastafilespath"};
		#tax=${tax%".best.ed.fasta"};
		#tax=${tax%".ed.fasta"};
		echo $taxon;
		exonerate --verbose 0 --model protein2genome $myref $taxon  --showvulgar no --showalignment no --ryo "$mygene,$taxon,%ql,%qal,%qab,%qae,%ti\n" >> ${cwd}/$mygene.results.csv;
		exonerate --verbose 0 --model protein2genome $myref $taxon  --showvulgar no --showalignment no --ryo ">$taxon,%ti,%qab\n%tcs\n" >> ${cwd}/$mygene.exons.fasta;
		LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d  $cwd\/$mygene.results.csv > $cwd\/$mygene.results.sorted.csv;
		#rm $mygene.results.csv;
	done	
fi

	### if the file matches the gene name then print it all out to an oputout fille
#	if [[$fastafilespath == *$mygene* ]] then

#### put it all in a loop two loops there will be... 

##
####step 3   ## create tests cases 
#echo "step 3: getting exons";
#echo " $cwd/bin/getcontigs.pl $cwd $myoverlap";
perl $cwd/bin/getcontigs.pl $cwd $myoverlap;
#echo "done with step 3";
##
####step 4
#echo "step 4, stichting";
echo $myoverlap $myaTRAMfile $cwd;
perl $cwd\/bin/stitch.contigs.pl  $myoverlap;
echo "done with step 4";
##
#####step 5
#echo "step 5, getting summary data";
#perl $cwd\/bin/summarystats.first.pl >> $summaryfile;
#echo "done with step 5";
##

#if $debug
if cmp -s GeneA.overlap.10.contig_list.csv test_suite/debug_contig_list ; then
	echo "Testing contig list:  passed test!"
else
	echo "Testing contig list: failed test need to go back and check code" 
fi
 
if cmp -s GeneA.stitched_exons.fasta test_suite/debug_stitched_contigs ; then
	echo "Testing contig stitching:  passed test!"
else
	echo "Testing contig stitching: failed test need to go back and check code" 
fi

###############################################################################
######### if not debug then remove the results.csv and the sorted and stats... 
##############################################################################


#exonerate --model protein2genome $myref $mygene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --ryo "$mygene,%ql,%qal,%qab,%qae,abyss,%ti\n" >> $mygene.exonerate2.out;
#$cwd\/bin/fix_csv.pl $mygene.exonerate2.out $mygene.exonerate2.csv; 
#LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $mygene.exonerate2.csv > $mygene.exonerate2.sorted.csv;
#exonerate --model protein2genome $myref $mygene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --verbose 0 --ryo ">$mygene,assembler,%ti,%qab,%qae\n%tcs\n" >> $mygene.exonerate2.fasta;
##
#perl $cwd/bin/editcontigs.pl $cwd;
#echo "done with step 6";
##
#####step 7
#echo "step 7, running second get contigs";
#perl $cwd\/bin/getcontigs.second.pl $cwd;
#echo "done with step 7";
##
#####step 8
#echo "step 8, stitching exonerate contigs";
#perl $cwd\/bin/stitch.Exonerate.contigs.pl;
#echo "done with step 8";
####
##
#####step 9
#echo "step 9, generating summary data for the final time";
#perl $cwd\/bin/summarystats.second.pl $cwd  >> $summaryfile;
#echo "done with step 9";
##
##
###### final cleanup
#rm *.exonerate*;
#rm *.results*
#rm *aTRAM.Exonerate.Round1.OVERLAP.10.fasta;
#rm files
#
#
############################################################################################################
