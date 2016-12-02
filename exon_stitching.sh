#!/usr/bin/bash

########################################################################################
#shell script to run post aTRAM pipeline created by J M Allen.
#script written by B M Boyd, UIUC
#Date created: 6-May-2014
_______________________________________________________________________________________
#USAGE:
#Usage: copy this shell script into your directory with best files
#Usage: use text editor  to input variables in block below
#Usage: sh Loop_Exonerate_pipeline.sh
_______________________________________________________________________________________


#VARIABLES TO DEFINE BEFORE RUNNING (substitute "<<repalce me>>" with your "variable"):
myaTRAMfile='out.best';  #### input common  file names for aTRAM contigs, eg. "I4.trinity.ed"

### Look through the whole thing and add the .ed
## THIS IS A BIT CONFUSING AS THE FILENAME DOES NOT EXIST UNTIL AFTER THE FIRST SCRIPT and should not include the .fasta.
summaryfile='SummaryFile.txt';  ### input a file name of your choice, this will create a table with statistics for all the taxa
#mypath='/srv/gspeedq32/Bestfiles_test_for_julie'; ### pathway to directory with the *.best.fasta files , leave off last "/"
pathtoreference='/data/kim/psocids/PHUM_mod_refs';  ##leave off the last /
myoverlap='10';  ### overlap to use, typically is "10" but can change if needed  #### THIS SHOULD BE 10 FOR THE FIRST ITERATION AND 0 FOR THE SECOND
########################################################################################

cwd=$(pwd);

#step 1
echo "step 1, edit best files";
perl $cwd\/bin/editbestfiles.pl $cwd;
echo "done with step 1";

#step 2
echo "step 2, running exonerate the first time";
for myref  in ${pathtoreference}/*.fasta; do 
	mygene=`expr match $myref '.*\(PHUM[0-9]*_EOG7[A-Z0-9]*\)'`;
	#echo $mygene;
############# NEED A LIST OF LIBRARIES ######
	for lib in  Anama Arnom; do   #       P01_WA1 P01_WA2 P01_WA3 P01_WA4 P01_WA5 P01_WA7 P01_WA8 P01_WA9 P01_WA10; do  
		exonerate --model protein2genome ${myref} ${cwd}/${mygene}.${lib}.${myaTRAMfile}.ed.fasta --showvulgar no --showalignment no --ryo "${mygene},${lib},%ql,%qal,%qab,%qae,abyss,%ti\n" >> ${cwd}/${mygene}.results.out;
	done;
#	###### need each lib to be done separately and printed to one file then run the fix.csv script.
	$cwd\/bin/fix_csv.pl $cwd\/$mygene.results.out $cwd\/$mygene.results.csv;
	LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $cwd\/$mygene.results.csv > $cwd\/$mygene.results.sorted.csv;
	done;
echo "done with step 2";

#
#step 3
echo "step 3, running first get contigs";
perl $cwd/bin/getcontigs.first.pl;
echo "done with step 3";
#
#step 4
echo "step 4, stichting";
echo $myoverlap $myaTRAMfile $cwd;
perl $cwd\/bin/stitch.aTRAM.contigs.pl  $myoverlap $myaTRAMfile $cwd;
echo "done with step 4";
#
##step 5
echo "step 5, getting summary data";
perl $cwd\/bin/summarystats.first.pl >> $summaryfile;
echo "done with step 5";
#
#step 6

#cd $pathtoreference
echo "step 6, running exonerate for second time";
for secref in ${pathtoreference}/*.fasta; do
	secgene=`expr match $secref '.*\(PHUM[0-9]*_EOG7[A-Z0-9]*\)'`;  
	echo $secref;
	exonerate --model protein2genome $secref $cwd/$secgene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --ryo "$secgene,%ql,%qal,%qab,%qae,abyss,%ti\n" >> $cwd/$secgene.exonerate2.out;
	$cwd\/bin/fix_csv.pl $cwd/$secgene.exonerate2.out $cwd/$secgene.exonerate2.csv; 
	LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $cwd/$secgene.exonerate2.csv > $cwd/$secgene.exonerate2.sorted.csv;
	exonerate --model protein2genome $secref $cwd/$secgene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --verbose 0 --ryo ">$secgene,abyss,%ti,%qab,%qae\n%tcs\n" >> $cwd/$secgene.exonerate2.fasta;
done;

perl $cwd/bin/editcontigs.pl $cwd;
echo "done with step 6";


##step 7
echo "step 7, running second get contigs";
perl $cwd\/bin/getcontigs.second.pl $cwd;
echo "done with step 7";
#
##step 8
echo "step 8, stitching exonerate contigs";
perl $cwd\/bin/stitch.Exonerate.contigs.pl;
echo "done with step 8";
#
##step 9
echo "step 9, generating summary data for the final time";
perl $cwd\/bin/summarystats.second.pl $cwd  >> $summaryfile;
echo "done with step 9";
