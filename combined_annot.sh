#!/bin/sh -l

strain_name=$1
out_dir=$2
fasta_file=$3
annot_ref=$4
annot_maker=$5
snap_dir=$6
SCRIPTS=$7
protein_file_type=$8 # SGD or ENSEMBL
cutoff=$9 
blast_db_dir=${10}

. $SCRIPTS/configs.cf

temp_dir=$out_dir/temp
rm -rf $temp_dir
mkdir -p $temp_dir

dir1=$annot_ref # annotations using lastz versus s288c
dir2=$annot_maker
codex=$dir1/$strain_name.codex
maker_gff=$dir2/genes.gff

blast_out=$out_dir/blast_out
gff=$out_dir/gff

mkdir -p $blast_out
mkdir -p $gff

rm -rf $out_dir/$strain_name.all.genes.gff
echo "#" > $out_dir/temp.codex

less $fasta_file | grep ">" | sed 's/>//g' > $out_dir/scf.list
while read scf_line 
do
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
	
	if [ $scf_len -gt 300 ]
	then
		rm -rf $temp_dir/temp*
		$BIN/extract_scaf_codex $codex $scf_name > $temp_dir/temp.codex
		less $maker_gff | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/temp.gff
		$BIN/merge_gff $temp_dir/temp.codex $temp_dir/temp.gff > $temp_dir/temp.genes.gff
		less $temp_dir/temp.genes.gff | grep -w 'gene\|match' | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/temp.gff
		$BIN/pull_fasta_scaf $fasta_file $scf_name > $temp_dir/temp.fasta
		$BIN/check_start_stop_codons $temp_dir/temp.fasta $temp_dir/temp.gff >> $temp_dir/temp.genes.gff
		$BIN/gff2codex $temp_dir/temp.genes.gff CDS > $temp_dir/temp.genes.codex
		$BIN/reverse_exon_order $temp_dir/temp.genes.codex > $temp_dir/temp1.genes.codex
		num_len=`less "$temp_dir"/temp1.genes.codex | wc -l | awk '{print $1}'`
		if [ $num_len -gt 0 ]
		then
			$BIN/pull_c $temp_dir/temp.fasta $temp_dir/temp1.genes.codex > $temp_dir/temp.genes.dna
			$BIN/dna2aa -v $temp_dir/temp.genes.dna 1 > $temp_dir/temp.genes.aa
			$BIN/check_aa $temp_dir/temp.genes.gff $temp_dir/temp.genes.aa > $temp_dir/temp.filtered.genes.gff
			$BIN/check_splice_signal $temp_dir/temp.fasta $temp_dir/temp.filtered.genes.gff > $temp_dir/temp1.filtered.genes.gff
			$BIN/merge_gff $out_dir/temp.codex $temp_dir/temp1.filtered.genes.gff > $temp_dir/temp.all.gff
			$BIN/gff2codex $temp_dir/temp.all.gff > $temp_dir/temp.all.codex
			num_len=`less "$temp_dir"/temp.all.codex | wc -l | awk '{print $1}'`
		fi
		if [ $num_len -gt 0 ]
		then
			$BIN/pull_c $temp_dir/temp.fasta $temp_dir/temp.all.codex > $temp_dir/temp.dna
      cd $blast_db_dir
      $BLAST/blastx -db ref -query $temp_dir/temp.dna -outfmt "7 sallacc pident evalue stitle" -out $temp_dir/temp.blastx.out
      cd $temp_dir
		  $BIN/update_gff_blastx $temp_dir/temp.all.gff $temp_dir/temp.blastx.out $protein_file_type $cutoff > $temp_dir/temp.all.genes.gff
			$BIN/rm_redun_gff $temp_dir/temp.all.genes.gff > $temp_dir/temp.final.genes.gff
			$BIN/merge_gff $out_dir/temp.codex $temp_dir/temp.final.genes.gff >> $gff/$strain_name.genes.gff
			less $temp_dir/temp.blastx.out >> $blast_out/$strain_name.blastx.out
		fi
	fi
done < $out_dir/scf.list
