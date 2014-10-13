#!/bin/sh -l

chainParams="-minScore=3000 -linearGap=loose"

axtChainNet=$1
seq_name=$2
out_dir=$3
ref_dir=$4
ref_name=$5
seq_file=$6
BIN=$7

#seq_file=$out_dir/$seq_name.fasta
main_dir=$out_dir/$ref_name.chain.net/$ref_name
target=$main_dir/$seq_name

mkdir -p $out_dir/$ref_name.chain.net
mkdir -p $out_dir/$ref_name.chain.net/$ref_name
mkdir -p $target
mkdir -p $target/lav
mkdir -p $target/chain
mkdir -p $target/net

if [ ! -e $main_dir/$seq_name.size ]
then
	$axtChainNet/faSize $seq_file -detailed > $main_dir/$seq_name.size
fi

for chr_seq in "$ref_dir"/*".fa" 
do
	chr_name=`basename $chr_seq | cut -d '.' -f1`
	echo $chr_name $chr_seq
	$BIN/lastz $chr_seq $seq_file --ambiguous=iupac --format=lav > $target/lav/$chr_name.$seq_name.lav
	$axtChainNet/lavToAxt $target/lav/$chr_name.$seq_name.lav -tfa $chr_seq -fa $seq_file $target/lav/$chr_name.$seq_name.axt
	$axtChainNet/axtChain $chainParams $target/lav/$chr_name.$seq_name.axt -faT $chr_seq -faQ $seq_file $target/chain/$chr_name.chain
#axtChainNet/faSize $ref1 -detailed > $main_dir/$chr_name.size
	$axtChainNet/chainNet $target/chain/$chr_name.chain -minSpace=1 $ref_dir/$chr_name.size $main_dir/$seq_name.size $target/net/$chr_name.temp.net $target/net/query.temp.net 
	$axtChainNet/netSyntenic $target/net/$chr_name.temp.net $target/net/$chr_name.net
done

rm -rf $target/net/*.temp.net
