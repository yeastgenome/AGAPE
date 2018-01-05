#!/bin/sh -l

fasta=$1
gff_file=$2
out_dir=$3
strain_name=$4
SCRIPTS=$5

. $SCRIPTS/configs.cf

mkdir -p $out_dir/temp
temp_dir=$out_dir/temp
mkdir -p $out_dir/gff
mkdir -p $out_dir/cds
mkdir -p $out_dir/aa

rm -rf $out_dir/gff/$strain_name.gff
rm -rf $out_dir/cds/$strain_name.cds.fasta
rm -rf $out_dir/aa/$strain_name.aa.fasta

echo "#" > $out_dir/temp.codex
less $fasta | grep ">" | sed 's/>//g'  > $out_dir/scf.list

count=1
while read scf_line
do
	rm -rf $temp_dir/cds.fasta
	rm -rf $temp_dir/scf.fasta
	rm -rf $temp_dir/temp.genes.gff
	rm -rf $temp_dir/temp.codex
	rm -rf $temp_dir/genes.codex
	rm -rf $temp_dir/genes.dna
	rm -rf $temp_dir/genes.aa
	rm -rf $temp_dir/temp.gff
	rm -rf $temp_dir/valid.genes.gff
	rm -rf $temp_dir/temp.ordered.codex
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
	if [ $scf_len -gt 300 ]
	then
		echo $scf_name
		less $gff_file | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/temp.genes.gff
		$BIN/gff2codex $temp_dir/temp.genes.gff CDS > $temp_dir/temp.codex
		$BIN/reverse_exon_order $temp_dir/temp.codex > $temp_dir/temp.ordered.codex
		$BIN/check_len_codex $temp_dir/temp.ordered.codex > $temp_dir/genes.codex
		num=`less $temp_dir/genes.codex | wc -l`
		if [ $num -gt 0 ]
		then
			$BIN/pull_fasta_scaf $fasta $scf_name > $temp_dir/scf.fasta
			$BIN/pull_c $temp_dir/scf.fasta $temp_dir/genes.codex > $temp_dir/genes.dna
			#BIN/merge_gff $temp_dir/genes.codex $out_dir/temp.gff > $temp_dir/temp.gff
			$BIN/dna2aa -v $temp_dir/genes.dna 1 > $temp_dir/genes.aa
			$BIN/check_aa $temp_dir/temp.genes.gff $temp_dir/genes.aa LAST_COLUMN | grep "gene\|CDS" > $temp_dir/valid.genes.gff
		fi
		num=`less $temp_dir/valid.genes.gff | grep -w "gene" | wc -l`
	else
		num=0
	fi

	if [ $num -gt 0 ]
	then
		num_genes=`less $temp_dir/valid.genes.gff | grep -w "gene" | wc -l`
		#   less $temp_dir/$scf_name.valid.genes.gff >> $out_dir/final_gff/$strain_name.non_ref.gff
		$BIN/gff2codex $temp_dir/valid.genes.gff CDS_NUM $count > $temp_dir/temp.codex
		less $temp_dir/valid.genes.gff >> $out_dir/gff/$strain_name.gff
		count=`expr $count + $num_genes`
		$BIN/pull_c $temp_dir/scf.fasta $temp_dir/temp.codex > $temp_dir/cds.fasta
		$BIN/dna2aa -v $temp_dir/cds.fasta 1 >> $out_dir/aa/$strain_name.aa.fasta
		less $temp_dir/cds.fasta >> $out_dir/cds/$strain_name.cds.fasta
	fi
done < $out_dir/scf.list

#mv $temp_dir/final.genes.gff $out_dir/$strain_name.gff
