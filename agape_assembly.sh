#!/bin/sh

# --- input file in FASTQ with the full path ---
# --- output file is written as $out_dir/$out_name.scf.fasta ---

SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
phred_type=33 # quality score type; change this to 64 for Illumina 1.3 and 1.5
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

#fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

mode=1 # 1: single end, 2: paired end

if [ $# -ne 4 ] && [ $# -ne 5 ]
then
  echo "Usage: $SCRIPT out_dir output_name AGAPE_main_path sequence1 [or sequence2 for paired end]"
  exit 1
fi

out_dir=$1
out_name=$2

SCRIPTS=$3

. $SCRIPTS/configs.cf

seq1=$4

temp_dir=$out_dir/"$out_name"_assembly
mkdir -p $temp_dir

if [ $# == 4 ]
then
  mode=1
	$SCRIPTS/error_correction.sh $out_dir $out_name $phred_type $seq1 $SCRIPTS # fastq files after error correction are named $out_dir/$out_name.*.fastq
	$SCRIPTS/assemble.sh $out_dir $out_name $mode $SCRIPTS
elif [ $# == 5 ]
then
  mode=2
  seq2=$5
	$SCRIPTS/error_correction.sh $out_dir $out_name $phred_type $seq1 $seq2 $SCRIPTS
	$SCRIPTS/assemble.sh $out_dir $out_name $mode $SCRIPTS
else
  echo "Usage: $SCRIPT output_name sequence1 [or sequence2 for paired end]"
  exit 1
fi
