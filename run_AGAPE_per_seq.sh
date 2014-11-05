#!/bin/sh

SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
BIN=/srv/gs1/projects/cherry/giltae/AGAPE/bin 
phred_type=33 # quality score type; change this to 64 for Illumina 1.3 and 1.5
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

#fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

. $SCRIPTS/configs.cf 

out_dir=$1
out_name=$2
seq1=$3

if [ $# == 3 ]
then
	$SCRIPTS/agape_assembly.sh $out_dir $out_name $SCRIPTS $seq1
	mode=1
elif [ $# == 4 ]
then
	seq2=$4
	mode=2
	$SCRIPTS/agape_assembly.sh $out_dir $out_name $SCRIPTS $seq1 $seq2
fi

contigs=$out_dir/$out_name.scf.fasta # assembly results from agape_assembly.sh
$SCRIPTS/agape_annot.sh $out_dir $out_name $contigs $SCRIPTS # resutls file in GFF is $out_dir/comb_annot/$out_name.gff

if [ $# == 3 ]
then
	$SCRIPTS/agape_novel_genes.sh $our_dir $out_name $SCRIPTS $out_dir/$out_name.scf.fasta $out_dir/comb_annot/$out_name.gff "$out_dir"/"$out_name"_1.pe.fastq 
elif [ $# == 4 ]
then
	$SCRIPTS/agape_novel_genes.sh $our_dir $out_name $SCRIPTS $out_dir/$out_name.scf.fasta $out_dir/comb_annot/$out_name.gff "$out_dir"/"$out_name"_1.pe.fastq "$out_dir"/"$out_name"_2.pe.fastq 
fi
# results file in GFF is $out_dir/non_ref/$out_name.novel.orfs.gff
