#!/bin/sh -l

#seq_dir=/srv/gs1/projects/cherry/giltae/pan_genome/fastq
seq_dir=/srv/gs1/projects/cherry/giltae/strains/validation/fastq
#out_dir=/srv/gs1/projects/cherry/giltae/strains/hugeseq
out_dir=/srv/gs1/projects/cherry/giltae/AGAPE/hugeseq
#list=/srv/gs1/projects/cherry/giltae/duke/hugeseq_rerun1.list
list=/srv/gs1/projects/cherry/giltae/strains/validation/fastq/liti_sc_list.txt

while read line
do
	file_name=`echo $line | awk '{print $1}'`
	cur_name=`echo $line | awk '{print $2}'`
	ln -s $seq_dir/"$file_name"_1.fastq $seq_dir/"$cur_name"_1.fastq
	ln -s $seq_dir/"$file_name"_2.fastq $seq_dir/"$cur_name"_2.fastq
	mkdir -p $out_dir/$cur_name
	cd $out_dir/$cur_name
	rm -rf $out_dir/$cur_name/*
	qsub -cwd -o $out_dir/$cur_name -e $out_dir/$cur_name -V -l h_vmem=12G -l h_stack=10M -q extended /srv/gs1/projects/cherry/giltae/AGAPE/hugeseq.sh $cur_name $seq_dir $out_dir/$cur_name
#	echo $cur_name $seq_dir $out_dir
done < $list
