#!/bin/sh

# --- input file in $out_dir/$out_name.scf.fasta ---

#SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
#SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

if [ $# -ne 4 ] && [ $# -ne 5 ]
then
  echo "Usage: agape_novel_genes.sh output_directory output_name AGAPE_main_path seq1 (or seq2)" 
  exit 1
fi

out_dir=$1
strain_name=$2
SCRIPTS=$3
seq1=$4

. $SCRIPTS/configs.cf

fasta=$out_dir/$strain_name.scf.fasta
gff=$out_dir/comb_annot/$strain_name.gff

non_ref_dir=$out_dir/non_ref
rm -rf $non_ref_dir
mkdir -p $non_ref_dir

if [ -f $REF_FASTA.ann ]
then
	echo "bwa index exists"
else
	$BWA/bwa index $REF_FASTA
fi

cd $non_ref_dir
if [ $# == 5 ]
then
#	seq2=$out_dir/"$strain_name"_2.pe.fastq
	seq2=$5
	$SCRIPTS/run_non_ref.sh $non_ref_dir $strain_name $SCRIPTS $fasta $gff $seq1 $seq2
elif [ $# == 4 ]
then
	$SCRIPTS/run_non_ref.sh $non_ref_dir $strain_name $SCRIPTS $fasta $gff $seq1
fi
