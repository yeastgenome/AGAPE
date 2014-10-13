#!/bin/sh -l

gff_file=$1
out_dir=$2
strain_name=$3
SCRIPTS=$4

mkdir -p $out_dir/temp
temp_dir=$out_dir/temp

rm -rf $out_dir/non_orf.fasta

echo "#" > $out_dir/temp.codex
less $fasta_file | grep ">" | sed 's/>//g'  > $out_dir/scf.list
while read scf_line
do
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
  if [ $scf_len -gt 300 ]
	then
		less $gff_file | grep -w 'gene' | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/temp.genes.gff
		$SCRIPTS/merge_gff $out_dir/temp.codex $temp_dir/temp.genes.gff >> $temp_dir/final.genes.gff
	fi
done < $out_dir/scf.list

mv $temp_dir/final.genes.gff $out_dir/$strain_name.gff
