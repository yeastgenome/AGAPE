#!/bin/sh

# --- input file in $out_dir/$out_name.scf.fasta ---

#SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

if [ $# -ne 6 ] && [ $# -ne 7 ]
then
  echo "Usage: $SCRIPT output_directory output_name AGAPE_main_path fasta gff sequence1 sequence2"
  exit 1
fi

out_dir=$1
strain_name=$2

SCRIPTS=$3

fasta=$4
gff=$5
seq1=$6

. $SCRIPTS/configs.cf

non_ref_dir=$out_dir/non_ref
mkdir -p $non_ref_dir

$BWA/bwa index $REF_FASTA
cd $non_ref_dir
if [ $# == 7 ]
then
	seq2=$7
	$SCRIPTS/run_non_ref.sh $non_ref_dir $strain_name $SCRIPTS $fasta $gff $seq1 $seq2
elif [ $# == 6 ]
then
	$SCRIPTS/run_non_ref.sh $non_ref_dir $strain_name $SCRIPTS $fasta $gff $seq1
fi
