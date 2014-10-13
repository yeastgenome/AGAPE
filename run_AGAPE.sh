#!/bin/sh

AGAPE_main_dir=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

for in `ls $fastq_dir`
do
	$AGAPE_main_dir/error_correction.sh 
done
