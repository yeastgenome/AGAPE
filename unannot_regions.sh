#!/bin/sh -l

fasta_file=$1
gff_file=$2
out_dir=$3
SCRIPTS=$4

. $SCRIPTS/configs.cf

mkdir -p $out_dir/temp_unannot

temp_dir=$out_dir/temp_unannot

rm -rf $out_dir/non_orf.fasta

echo "#" > $out_dir/temp.codex
less $fasta_file | grep ">" | sed 's/>//g'  > $out_dir/scf.list
while read scf_line
do
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
  if [ $scf_len -gt 300 ]
	then
		$BIN/pull_fasta_scaf $fasta_file $scf_name > $temp_dir/temp.fasta
		head -1 $temp_dir/temp.fasta > $temp_dir/head.txt
		less $gff_file | grep -w 'gene' | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/temp.genes.gff
		$BIN/merge_gff $out_dir/temp.codex $temp_dir/temp.genes.gff > $temp_dir/known.genes.gff
		$BIN/non_ref $temp_dir/head.txt $temp_dir/known.genes.gff > $temp_dir/non_orf.codex

		num_len=`less "$temp_dir"/non_orf.codex | wc -l | awk '{print $1}'`

		if [ $num_len -gt 0 ]
		then
			$BIN/pull_c $temp_dir/temp.fasta $temp_dir/non_orf.codex >> $out_dir/non_orf.fasta
		fi
	fi
done < $out_dir/scf.list

rm -rf $temp_dir
