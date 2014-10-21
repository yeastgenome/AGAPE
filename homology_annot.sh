#!/bin/sh -l

SCRIPTS=$3

if [ ! -e $SCRIPTS/configs.cf ]
then
	echo "$SCRIPTS/configs.cf not exist"
	exit 1
fi

. $SCRIPTS/configs.cf

CHR_SEQ_DIR=$REF_DIR/$REF_NAME.seq

if [ ! -d "$REF_DIR" ]
then
	echo "$REF_DIR not defined or in wrong path in ../configs.cf"
	exit 1
fi

if [ -z "$REF_NAME" ]
then
	echo "REF_NAME not assigned in ../configs.cf"
	exit 1
fi

cur_dir=$1
seq_name=$2
annot_dir=$cur_dir/annot

mkdir -p $CHR_SEQ_DIR
mkdir -p $annot_dir

if [ ! -e $REF_DIR/$REF_NAME.size ]
then
	if [ ! -e $REF_FASTA ]
	then
		echo "error: $REF_FASTA not exist or moved to other location"
	fi
	$axtChainNet/faSize $REF_FASTA -detailed > $REF_DIR/$REF_NAME.size
fi

mkdir -p $REF_DIR/chr_seq

for chr_name in `less "$REF_DIR/$REF_NAME.size" | awk '{print $1}'`
do
	$BIN/pull_fasta_scaf $REF_FASTA $chr_name > $REF_DIR/chr_seq/$chr_name.fa
	$axtChainNet/faSize $REF_DIR/chr_seq/$chr_name.fa -detailed > $REF_DIR/chr_seq/$chr_name.size
done

if [ ! -e "$cur_dir"/"$seq_name".scf.fasta ]
then
	echo "$cur_dir/$seq_name.scf.fasta moved or not exist"
	exit 1
fi

ln -s $cur_dir/$seq_name.scf.fasta $annot_dir/$seq_name.fasta
$SCRIPTS/ChainNet.sh $axtChainNet $seq_name $annot_dir $REF_DIR/chr_seq $REF_NAME $annot_dir/$seq_name.fasta $BIN #three result directories lav, chain, and net in $annot_dir/$REF_NAME.chain.net
$SCRIPTS/intervals.sh $seq_name $annot_dir $REF_DIR $REF_NAME $annot_dir/$REF_NAME.chain.net $BIN #$REF_DIR/intervals are created

rm -rf $annot_dir/$seq_name.codex
while read line
do
	chr_name=`echo $line | awk '{print $1}'`
	rm -rf $annot_dir/codex/$chr_name.codex
	$SCRIPTS/infer-scaf-annot.sh $annot_dir $REF_DIR $REF_NAME $chr_name $BIN $SCRIPTS
	cat $annot_dir/codex/$chr_name.codex >> $annot_dir/$seq_name.codex
done < $REF_DIR/$REF_NAME.size
