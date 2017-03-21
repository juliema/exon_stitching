#!/usr/bin/sh

########################################################################################
#Usage: sh exon_stitching.sh AminoAcidFile.fasta path_to_target_assemblies/ output_name
########################################################################################

# help menu
if [ "$1" = "-h" ] ; then
  echo "Usage: `basename $0` AminoAcidFile.fasta  path_to_target_assemblies/ ";
  exit 0
fi

echo amino acid file $1;
echo path to atram assemblies $2; 
myref=$1;
bestfilespath=$2;
summaryfile='SummaryFile.txt';  ### input a file name of your choice, this will create a table with statistics for all the taxa
myoverlap='10';  ### overlap to use, typically is "10" but can change if needed  #### THIS SHOULD BE 10 FOR THE FIRST ITERATION AND 0 FOR THE SECOND
cwd=$(pwd);
mygene=${myref%.fasta} 
#mygene=`echo $myref | egrep -Eo  '^[^.fasta]*'`;
echo $mygene;

#step 1
echo "step 1, edit best files";
echo targeting the following atram assemblies;
perl $cwd\/bin/editbestfiles.pl $bestfilespath;
echo "done with step 1";

#step 2
echo "step 2, running exonerate the first time";
for lib in $bestfilespath*.ed.fasta; do
	echo $lib;
	exonerate --model protein2genome ${myref} $lib  --showvulgar no --showalignment no --ryo "${myref},${lib},%ql,%qal,%qab,%qae,abyss,%ti\n" >> ${cwd}/$mygene.results.out;
	$cwd\/bin/fix_csv.pl $cwd\/$mygene.results.out $cwd\/$mygene.results.csv;
	LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $cwd\/$mygene.results.csv > $cwd\/$mygene.results.sorted.csv;
done;
#
###step 3   ## create tests cases 
echo "step 3, running first get contigs";
perl $cwd/bin/getcontigs.first.pl;
echo "done with step 3";
#
###step 4
echo "step 4, stichting";
echo $myoverlap $myaTRAMfile $cwd;
perl $cwd\/bin/stitch.aTRAM.contigs.pl  $myoverlap $myaTRAMfile $cwd;
echo "done with step 4";
#
####step 5
echo "step 5, getting summary data";
perl $cwd\/bin/summarystats.first.pl >> $summaryfile;
echo "done with step 5";
#
###step 6
exonerate --model protein2genome $myref $mygene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --ryo "$mygene,%ql,%qal,%qab,%qae,abyss,%ti\n" >> $mygene.exonerate2.out;
$cwd\/bin/fix_csv.pl $mygene.exonerate2.out $mygene.exonerate2.csv; 
LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $mygene.exonerate2.csv > $mygene.exonerate2.sorted.csv;
exonerate --model protein2genome $myref $mygene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --verbose 0 --ryo ">$mygene,assembler,%ti,%qab,%qae\n%tcs\n" >> $mygene.exonerate2.fasta;
#
perl $cwd/bin/editcontigs.pl $cwd;
echo "done with step 6";
#
####step 7
echo "step 7, running second get contigs";
perl $cwd\/bin/getcontigs.second.pl $cwd;
echo "done with step 7";
#
####step 8
echo "step 8, stitching exonerate contigs";
perl $cwd\/bin/stitch.Exonerate.contigs.pl;
echo "done with step 8";
###
#
####step 9
echo "step 9, generating summary data for the final time";
perl $cwd\/bin/summarystats.second.pl $cwd  >> $summaryfile;
echo "done with step 9";
#
#
##### final cleanup
rm *.exonerate*;
rm *.results*
rm *aTRAM.Exonerate.Round1.OVERLAP.10.fasta;
rm files



