#!/bin/sh

SCRIPTS=$1
. $SCRIPTS/configs.cf

TEMP=$REF_DIR/$REF_NAME.temp
mkdir -p $TEMP

while read line
do
  chr_name=`echo $line | awk '{print $1}'`
	$BIN/gff2temp $REF_DIR/$REF_NAME.gff $chr_name > $TEMP/$chr_name.m_temp
	$BIN/get_gene_bound $TEMP/$chr_name.m_temp > $TEMP/$chr_name.temp_loc
	$BIN/pull_c $REF_DIR/chr_seq/$chr_name.fa $TEMP/$chr_name.temp_loc > $TEMP/$chr_name.temp_gene    # first nt of seq counts as "1"
	$BIN/pull_c $REF_DIR/chr_seq/$chr_name.fa $TEMP/$chr_name.m_temp > $TEMP/$chr_name.temp_dna    # first nt of seq counts as "1"
	$BIN/dna2aa -v $TEMP/$chr_name.temp_dna 1 > $TEMP/$chr_name.temp_aa
done < $REF_DIR/$REF_NAME.size
