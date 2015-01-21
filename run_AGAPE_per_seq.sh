#!/bin/sh

#BIN=/srv/gs1/projects/cherry/giltae/AGAPE/bin 
phred_type=33 # quality score type; change this to 64 for Illumina 1.3 and 1.5

#SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
#fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

out_dir=$1
out_name=$2
SCRIPTS=$3
seq1=$4

. $SCRIPTS/configs.cf 

if [ $# == 4 ]
then
	mode=1
	$SCRIPTS/agape_assembly.sh $out_dir $out_name $SCRIPTS $seq1
elif [ $# == 5 ]
then
	seq2=$5
	mode=2
	$SCRIPTS/agape_assembly.sh $out_dir $out_name $SCRIPTS $seq1 $seq2
fi

contigs=$out_dir/$out_name.scf.fasta # assembly results from agape_assembly.sh
$SCRIPTS/agape_annot.sh $out_dir $out_name $contigs $SCRIPTS # resutls file in GFF is $out_dir/comb_annot/$out_name.gff

if [ $# == 4 ]
then
	$SCRIPTS/agape_novel_genes.sh $out_dir $out_name $SCRIPTS $seq1
#out_dir"/"$out_name".scf.fasta $out_dir/comb_annot/"$out_name".gff "$out_dir"/"$out_name"_1.pe.fastq 
elif [ $# == 5 ]
then
	$SCRIPTS/agape_novel_genes.sh $out_dir $out_name $SCRIPTS $seq1 $seq2
#out_dir"/"$out_name".scf.fasta "$out_dir"/comb_annot/"$out_name".gff "$out_dir"/"$out_name"_1.pe.fastq "$out_dir"/"$out_name"_2.pe.fastq 
fi
# results file in GFF is $out_dir/non_ref/$out_name.novel.orfs.gff
