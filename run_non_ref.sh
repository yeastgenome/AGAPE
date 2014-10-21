#!/bin/sh

# --- input file in $out_dir/$out_name.scf.fasta ---

#SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
#BIN=/srv/gs1/projects/cherry/giltae/AGAPE/bin
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

if [ $# -ne 6 ] && [ $# -ne 7 ]
then
  echo "Usage: $SCRIPT output_directory output_name AGAPE_main_path assembly_contigs sequence1 sequence2"
  exit 1
fi

out_dir=$1
seq_name=$2

SCRIPTS=$3
. $SCRIPTS/configs.cf

assem_fasta=$4
gff=$5
seq1=$6

cd $out_dir
if [ $# == 7 ]
then
	seq2=$7
	$SCRIPTS/non_ref.sh $out_dir $seq_name $SCRIPTS $seq1 $seq2
	$SCRIPTS/assemble.sh $out_dir $seq_name 4 $SCRIPTS  
else
	$SCRIPTS/non_ref.sh $out_dir $seq_name $SCRIPTS $seq1
	$SCRIPTS/assemble.sh $out_dir $seq_name 3 $SCRIPTS
fi

## assembly results are in $out_dir/$seq_name.scf.fasta

$SCRIPTS/non_ref_contigs.sh $out_dir $seq_name $SCRIPTS $out_dir/$seq_name.scf.fasta $assem_fasta # resutls in $out_dir/$seq_name.final.inserted.assembly.intervals
$SCRIPTS/novel_orfs.sh $out_dir $seq_name $SCRIPTS $gff $out_dir/$seq_name.final.inserted.original.intervals $assem_fasta # results in $out_dir/$seq_name.novel.orfs.gff
