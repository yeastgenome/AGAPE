#!/bin/sh

<<<<<<< HEAD
AGAPE_main_dir=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

for in `ls $fastq_dir`
do
	$AGAPE_main_dir/error_correction.sh 
done
=======
SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path

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

	qsub -cwd -o $out_dir -e $out_dir -V -l h_vmem=6G -l h_stack=10M -q extended $SCRIPTS/run_AGAPE_per_seq.sh $out_dir $out_name $seq1 $seq2
done < /srv/gs1/projects/cherry/giltae/duke/strains.list
# results file in GFF is $out_dir/non_ref/$out_name.novel.orfs.gff
>>>>>>> 209b092fce7cd388a8c2e53b71d989a7bddd01b6
