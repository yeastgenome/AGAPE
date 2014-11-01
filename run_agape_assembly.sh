#!/bin/sh

SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
BIN=/srv/gs1/projects/cherry/giltae/AGAPE/bin 
phred_type=33 # quality score type; change this to 64 for Illumina 1.3 and 1.5
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

#fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

. $SCRIPTS/configs.cf 

seq_dir=/srv/gs1/projects/cherry/giltae/pan_genome/fastq

while read line
do
	out_name=`echo $line | awk '{print $1}'`
	out_dir=$SCRIPTS/duke/$out_name
	mkdir -p $out_dir
	cd $out_dir
	seq1=$seq_dir/"$out_name"_1.fastq
	seq2=$seq_dir/"$out_name"_2.fastq

	qsub -cwd -o $out_dir -e $out_dir -V -l h_vmem=6G -l h_stack=10M -q extended $SCRIPTS/agape_assembly.sh $out_dir $out_name $SCRIPTS $seq1 $seq2
done < /srv/gs1/projects/cherry/giltae/duke/rerun.list
# results file in GFF is $out_dir/non_ref/$out_name.novel.orfs.gff
