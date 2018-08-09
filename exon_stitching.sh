#!/bin/bash

########################################################################################
#Usage: sh exon_stitching.sh AminoAcidFile.fasta path_to_target_assemblies/ output_name
########################################################################################

# help menu
if [ "$1" = "-h" ] ; then
echo "Usage: bash `basename $0` list_of_taxa.txt reference_amino_acid_sequences.fasta  path_to_target_assemblies/  overlap";
exit 0
fi

#echo taxon list $1;
#echo amino acid reference file $2;
#echo path to assemblies $3;
#echo overlap $4;


######  Uncomment when not testing
mytax=$1;
myref=$2;
fastafilespath=$3;
#myoverlap=$4;  
cwd=$(pwd);
myoverlap=10;  


echo taxon list $mytax;
echo amino acid reference file $myref;
echo path to assemblies $fastafilespath;
echo overlap $myoverlap;

### First we edit all of the atram assemblies, remove empty files, get rid of 
##   weird haracters in the names and add a unique number to the end of each contig.
##   creates a file called Gene_list.txt with the genes and greates one fasta file for each reference gene.
perl $cwd\/bin/editfasta.pl $fastafilespath $myref;

### reads through the gene list and fines all the assembly files with that gene name in it and makes a list of all the assembly files for each gene
gene_file="Gene_list.txt";
while read -r line  || [[ -n "$line" ]]; do
	gene=$line
	gene=${gene%".reference.fasta"};
	echo "each gene is - $gene"
	find $fastafilespath -name "*$gene*.ed.fasta" >$gene.list.txt
done < "$gene_file"

#In this part we read through each reference then each taxon to run exonerate.
# we create a results file with the exonerate information for each exon
# we sort the results file by taxon and by location of the exon
# we use exonerate to pull out the exons
while read -r reffile || [[ -n "$reffile" ]]; do
	gene=${reffile%".reference.fasta"};
	echo ""
	echo reference file  $reffile;
	#echo $gene;
	while read -r taxfile || [[ -n "$taxfile" ]]; do
		tax=$taxfile
		echo taxon $tax
		while read -r assembly || [[ -n "$assembly" ]]; do
			if [[ $assembly  = *$tax.* ]]; then
				echo $assembly contains $tax
				exonerate --verbose 0 --model protein2genome $reffile $assembly  --showvulgar no --showalignment no --ryo "$gene,$tax,%ql,%qal,%qab,%qae,%ti\n" >> ${cwd}/$gene.results.csv;
				exonerate --verbose 0 --model protein2genome $reffile $assembly  --showvulgar no --showalignment no --ryo ">$tax,%ti,%qab\n%tcs\n" >> ${cwd}/$gene.exons.fasta;
				LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d  $cwd\/$gene.results.csv > $cwd\/$gene.results.sorted.csv;
                   		##rm $ref.results.csv;
			fi
		done < "$gene.list.txt"
	done < "$mytax"
done < "Gene_list.txt"


perl $cwd/bin/getcontigs.pl $cwd $myoverlap;
perl $cwd\/bin/stitch.contigs.pl  $myoverlap;
perl $cwd\/bin/summary_stats.pl $cwd;

#rm contig_files.txt;
#rm Gene_list.txt;
#rm *.sorted.csv;
#rm files;
#rm *results.csv
#rm $fastafilespath/*.ed.fasta;
#rm *.reference.fasta;

#### new comparison
 #  cmd  gene1.stitched_exons.fasta  test_suite/debug_stitched_contigs.fasta 
 #  diff  gene1.stitched_exons.fasta  test_suite/debug_stitched_contigs 
 #  cmd  gene1.overlap.10.contig_list.csv  test_suite/debug_gene1_contig_list.txt
 #  diff  gene1.overlap.10.contig_list.csv  test_suite/debug_gene1_contig_list.txt


#if $debug
#	echo "Testing one gene";
#  run as   bash exon_stitching.sh  GeneA.fasta  test_suite/ 10
#if cmp -s PHUM336790-PA.overlap.10.contig_list.csv test_suite/debug_contig_list ; then
#	echo "Testing contig determination:  passed test!"
#else
#	echo "Testing contig determination: failed test need to go back and check the get_contig section" 
#fi
#	if cmp -s PHUM336790-PA.stitched_exons.fasta test_suite/debug_stitched_contigs ; then
#		echo "Testing contig stitching:  passed test!";
#	else
#		echo "Testing contig stitching: failed test need to go back and check stitching contigs.";
#	fi
#fi

