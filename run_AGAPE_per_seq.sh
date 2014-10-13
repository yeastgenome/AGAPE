#!/bin/sh

SCRIPTS=/srv/gs1/projects/cherry/giltae/AGAPE # AGAPE main directory path
BIN=/srv/gs1/projects/cherry/giltae/AGAPE/bin 
out_dir=.
phred_type=33 # quality score type; change this to 64 for Illumina 1.3 and 1.5
SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

#fastq_dir=/srv/gs1/projects/cherry/giltae/AGAPE/output/fastq

. $SCRIPTS/configs.cf 

mode=1 # 1: single end, 2: paired end

if [ $# -ne 3 ] && [ $# -ne 4 ]
then
	echo "Usage: $SCRIPT output_name sequence1 [or sequence2 for paired end]"
	exit 1
fi 

out_name=$1
seq1=$2

if [ $# == 3 ] 
then
	mode=1
#SCRIPTS/error_correction.sh $out_dir $out_name $phred_type $seq1 $SCRIPTS
#SCRIPTS/assemble.sh $out_dir $out_name $mode $SCRIPTS
elif [ $# == 4 ]
then
	mode=2
	seq2=$3
#SCRIPTS/error_correction.sh $out_dir $out_name $phred_type $seq1 $seq2 $SCRIPTS
#SCRIPTS/assemble.sh $out_dir $out_name $mode $SCRIPTS
else 
	echo "Usage: $SCRIPT output_name sequence1 [or sequence2 for paired end]"
	exit 1
fi

if

if [ ! -e "$REF_FASTA" ]
then
	echo "No reference sequence data: annotation based on homology skipped"
else
	#SCRIPTS/homology_annot.sh $out_dir $out_name $SCRIPTS # results place in $out_dir/annot/$out_name.codex
fi

snap_dir=$out_dir/snap_files
mkdir -p $snap_files 

ln -s $REF_FASTA $snap_dir/$REF_NAME.dna
#SCRIPTS/prep_maker.sh $snap_dir $SCRIPTS

maker_dir=$out_dir/maker
rm -rf $out_dir/maker
mkdir -p $maker_dir

#SCRIPTS/run_maker.sh $out_dir $out_name $maker_dir $snap_dir $SCRIPTS # results in $maker_dir/genes.gff

comb_annot=$out_dir/comb_annot
rm -rf $comb_annot
mkdir -p $comb_annot

cd $comb_annot
$SCRIPTS/run_comb_annot.sh $comb_annot $BLAST $out_name $maker_dir $out_dir $snap_dir $SCRIPTS 

$SCRIPTS/final_annot.sh $comb_annot/gff/$out_name.genes.gff $comb_annot $out_name $SCRIPTS # resutls in $comb_annot/non_orf.fasta
