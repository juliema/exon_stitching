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


perl $cwd\/bin/editfasta.pl $fastafilespath $myref;

##step 2
file="Gene_list.txt";
cat $file | while read myref; do
	ref=${myref%".fasta"};
	echo $ref;
	for taxon in $fastafilespath*.ed.fasta; do
		echo $taxon;
		tax=${tax%".ed.fasta"};
		#echo $tax;
		exonerate --verbose 0 --model protein2genome $myref $taxon  --showvulgar no --showalignment no --ryo "$ref,$taxon,%ql,%qal,%qab,%qae,%ti\n" >> ${cwd}/$ref.results.csv;
		exonerate --verbose 0 --model protein2genome $myref $taxon  --showvulgar no --showalignment no --ryo ">$taxon,%ti,%qab\n%tcs\n" >> ${cwd}/$ref.exons.fasta;
		LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d  $cwd\/$ref.results.csv > $cwd\/$ref.results.sorted.csv;
		#rm $ref.results.csv;
	done
done

perl $cwd/bin/getcontigs.pl $cwd $myoverlap;
perl $cwd\/bin/stitch.contigs.pl  $myoverlap;
perl $cwd\/bin/summary_stats.pl Summary_stats.per.gene.csv;

rm contig_files.txt;
rm Gene_list.txt;
rm *.sorted.csv;
rm files;
rm *results.csv
rm test_suite/*.ed.fasta



#if $debug

echo "Testing one gene";
#  run as   bash exon_stitching.sh  GeneA.fasta  test_suite/ 10
#if cmp -s PHUM336790-PA.overlap.10.contig_list.csv test_suite/debug_contig_list ; then
#	echo "Testing contig determination:  passed test!"
#else
#	echo "Testing contig determination: failed test need to go back and check the get_contig section" 
#fi
 
if cmp -s PHUM336790-PA.stitched_exons.fasta test_suite/debug_stitched_contigs ; then
	echo "Testing contig stitching:  passed test!"
else
	echo "Testing contig stitching: failed test need to go back and check stitching contigs." 
fi

echo "Testing multiple genes";
# run as bash exon_stitching.sh Mult_genes.fasta test_suite/ 10



