#!/bin/sh

out_dir=$1
seq_name=$2

SCRIPTS=$3
. $SCRIPTS/configs.cf

gff=$4 # annotation in GFF 
interval_file=$5 # non-reference intervals 
fasta=$6 # assembly in FASTA

if [ -f $interval_file ]
then
  num_len=`less $interval_file | wc -l`
  if [ $num_len -gt 0 ]
  then
    while read line
    do
      scf_name=`echo $line | awk '{print $1}'`
      b=`echo $line | awk '{print $2}'`
      e=`echo $line | awk '{print $3}'`
      $BIN/gff_subset $gff $scf_name $b $e >> $out_dir/$seq_name.novel.orfs.gff
    done < $interval_file

    if [ -f $out_dir/$seq_name.novel.orfs.gff ]
    then
      num_len=`less $out_dir/$seq_name.novel.orfs.gff | wc -l`
      if [ $num_len -gt 0 ]
      then
        less $out_dir/$seq_name.novel.orfs.gff | grep -w "gene" > $out_dir/temp.gff
        count=1
        while read line
        do
          scf_name=`echo $line | awk '{print $1}'`
          b=`echo $line | awk '{print $4}'`
          e=`echo $line | awk '{print $5}'`
          $BIN/pull_fasta_scaf $fasta $scf_name > $out_dir/temp.fasta
          echo ">$seq_name.ORF$count $scf_name:$b-$e" >> $out_dir/$seq_name.novel.orfs.fasta
          $BIN/dna $b,$e $out_dir/temp.fasta | tail -n +2 >> $out_dir/$seq_name.novel.orfs.fasta
          count=`expr $count + 1`
        done < $out_dir/temp.gff
#       less $out_dir/$seq_name.novel.orfs.fasta >> $novel_orf_dir/novel.orfs.fasta
      fi
    fi
  fi
fi
