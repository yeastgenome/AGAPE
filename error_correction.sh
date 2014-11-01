#!/bin/sh -l


CORRECTION_K=41
out_dir=$1
seq_name=$2
phred_option=$3
seq1=$4
mode=1 # single-end, 2 for paired-end

SCRIPT=`echo $0 | sed -e 's;.*/;;'` # script name from command line; path removed for msgs

if [ $# == 5 ]
then
	mode=1
	SCRIPTS=$5
elif [ $# == 6 ]
then 
	mode=2
	seq2=$5
	SCRIPTS=$6
else # if [ $# -ne 5 ] && [ $# -ne 6 ]
  echo "Usage: $SCRIPT output_dir output_name phred_type sequence1 [or sequence2 for paired end]"
  exit 1
fi

echo "num of param: $#"

. $SCRIPTS/configs.cf
sga=$SGA/sga

rm -rf $out_dir/$seq_name.pp.ec.fa
rm -rf "$out_dir/$seq_name"_1.pe.fastq
rm -rf "$out_dir/$seq_name"_2.pe.fastq
rm -rf "$out_dir/$seq_name"_1.se.fastq
rm -rf "$out_dir/$seq_name"_2.se.fastq

fname1=$out_dir/"$seq_name"_1.sra.fastq
fname2=$out_dir/"$seq_name"_2.sra.fastq

if [ -e "$fname1" ]
then
	echo "$fname1 exists; change $fname1 to a different name"
	exit 1
fi

if [ ! -e "$SGA/sga" ]
then
	echo "$SGA/sga not installed or path not set"
	exit 1
fi

version=`less $seq1 | head -1 | awk '{print NF}'`

if [ "$version" == 1 ]
then
	if [ ! -f $seq1 ]
	then
		echo "$seq1 does not exist"
		exit 1
	else		
		ln -s $seq1 $fname1
	fi

	if [ $mode == 2 ]
	then
		if [ ! -f $seq2 ]
		then
			echo "$seq2 does not exist"
		else 
			ln -s $seq2 $fname2
		fi
	fi
else
	$BIN/edit_fq_head $seq1 1 > $fname1 # SGA and ABYSS require that header
	if [ $mode == 2 ]
	then
		$BIN/edit_fq_head $seq2 2 > $fname2
	fi
fi

#ln -s $seqdir/"$seq_name"_1.fastq $fname1
#ln -s $seqdir/"$seq_name"_2.fastq $fname2

if [ -e "$out_dir/$seq_name.pp.fa" ]
then
	echo "$out_dir/$seq_name.pp.fa exists"
elif [ "$phred_option" == 33 ]
then
	if [ $mode == 1 ]
	then
		$sga preprocess $fname1 > $out_dir/$seq_name.pp.fa
	else
		$sga preprocess -p 1 $fname1 $fname2 > $out_dir/$seq_name.pp.fa
	fi
else
	if [ $mode == 1 ]
	then
		$sga preprocess --phred64 $fname1 > $out_dir/$seq_name.pp.fa
	else
		$sga preprocess -p 1 --phred64 $fname1 $fname2 > $out_dir/$seq_name.pp.fa
	fi
fi

$sga index -a ropebwt --no-reverse $out_dir/$seq_name.pp.fa
$sga correct -k $CORRECTION_K --discard --learn $out_dir/$seq_name.pp.fa -o $out_dir/$seq_name.pp.ec.fa
#sga index -a ropebwt -t $CPU --no-reverse $out_dir/$seq_name.pp.fa
#sga correct -k $CORRECTION_K --discard --learn -t $CPU $out_dir/$seq_name.pp.fa -o $out_dir/$seq_name.pp.ec.fa

if [ $mode == 1 ]
then
	ln -s $out_dir/$seq_name.pp.ec.fa $out_dir/"$seq_name"_1.se.fastq
elif [ $version == 1 ]
then
	$BIN/separate_fq $out_dir/$seq_name.pp.ec.fa $out_dir/"$seq_name"_1.pe.fastq $out_dir/"$seq_name"_2.pe.fastq $out_dir/"$seq_name"_1.se.fastq $out_dir/"$seq_name"_2.se.fastq VERSION1
else
	$BIN/separate_fq $out_dir/$seq_name.pp.ec.fa $out_dir/"$seq_name"_1.pe.fastq $out_dir/"$seq_name"_2.pe.fastq $out_dir/"$seq_name"_1.se.fastq $out_dir/"$seq_name"_2.se.fastq 
fi

rm -rf $fname1
rm -rf $fname2
rm -rf $out_dir/$seq_name.pp.bwt
rm -rf $out_dir/$seq_name.pp.sai
rm -rf $out_dir/$seq_name.pp.fa
rm -rf $out_dir/$seq_name.pp.discard.fa
