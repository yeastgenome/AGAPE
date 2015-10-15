#!/bin/sh

ratio=0.5

out_dir=$1
seq_name=$2

SCRIPTS=$3
. $SCRIPTS/configs.cf

fasta=$4 # contigs of unmapped reads
assem_fasta=$5 # whole genome sequence assembly 

rm -rf $out_dir/$seq_name.inserted.assembly.intervals

for scf_name in `less "$assem_fasta" | grep ">" | awk '{print $1}' | cut -d '>' -f2`
do
  $BIN/pull_fasta_scaf $assem_fasta $scf_name > $out_dir/temp.fasta
  len=`$BIN/seq_len $out_dir/temp.fasta`
  lastz T=2 Y=3400 $out_dir/temp.fasta $REF_FASTA --ambiguous=iupac --format=maf > $out_dir/temp.maf
  $BIN/find_inserted_intervals $out_dir/temp.maf $scf_name $len >> $out_dir/$seq_name.inserted.assembly.intervals
  rm -rf $out_dir/temp.fasta
  rm -rf $out_dir/temp.maf
done

rm -rf $out_dir/$seq_name.inserted.intervals

for scf_name in `less "$fasta" | grep ">" | awk '{print $1}' | cut -d '>' -f2`
do
 $BIN/pull_fasta_scaf $fasta $scf_name > $out_dir/temp.fasta
 len=`$BIN/seq_len $out_dir/temp.fasta`
 $BIN/lastz T=2 Y=3400 $out_dir/temp.fasta $REF_FASTA --ambiguous=iupac --format=maf > $out_dir/temp.maf
 $BIN/find_inserted_intervals $out_dir/temp.maf $scf_name $len >> $out_dir/$seq_name.inserted.intervals
 rm -rf $out_dir/temp.maf
done

rm -rf $out_dir/temp.*
rm -rf $out_dir/$seq_name.inserted.original.intervals
rm -rf $out_dir/$seq_name.inserted.original.temp.intervals

while read line
do
 scf_name=`echo $line | awk '{print $1}'`
 b=`echo $line | awk '{print $2}'`
 e=`echo $line | awk '{print $3}'`

 $BIN/pull_fasta_scaf $fasta $scf_name > $out_dir/temp.seq.fasta
 $BIN/dna $b,$e $out_dir/temp.seq.fasta > $out_dir/temp.cur.seq.fasta
 $BIN/lastz T=2 Y=3400 $assem_fasta[multi] $out_dir/temp.cur.seq.fasta --ambiguous=iupac --format=maf > $out_dir/temp.seq.maf
 $BIN/scf_lift_over $out_dir/temp.seq.maf >> $out_dir/$seq_name.inserted.original.intervals
 rm -rf $out_dir/temp.seq.maf
 rm -rf $out_dir/temp.seq.fasta
 rm -rf $out_dir/temp.cur.seq.fasta
done < $out_dir/$seq_name.inserted.intervals

if [ -f $out_dir/$seq_name.inserted.original.intervals ]
then
 sort $out_dir/$seq_name.inserted.original.intervals > $out_dir/$seq_name.inserted.original.temp.intervals
 $BIN/merge_scf_intervals $out_dir/$seq_name.inserted.original.temp.intervals > $out_dir/$seq_name.inserted.original.intervals
else
 echo -n "" > $out_dir/$seq_name.inserted.original.intervals
fi

rm -rf $out_dir/$seq_name.final.inserted.original.intervals
if [ -f $out_dir/$seq_name.inserted.original.intervals ]
then
  num_len=`less $out_dir/$seq_name.inserted.original.intervals | wc -l`
  if [ $num_len -gt 0 ]
  then
    if [ -f $out_dir/$seq_name.inserted.assembly.intervals ]
    then
      num_len=`less $out_dir/$seq_name.inserted.assembly.intervals | wc -l`
      if [ $num_len -gt 0 ]
      then
        $BIN/common_scf_intervals $out_dir/$seq_name.inserted.original.intervals $out_dir/$seq_name.inserted.assembly.intervals > $out_dir/$seq_name.final.inserted.original.intervals
      fi
    fi
  fi
fi

