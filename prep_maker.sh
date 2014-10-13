#!/bin/sh -l

SCRIPT=`echo $0 | sed -e 's;.*/;;'`

if [ $# -ne 2 ]
then
	echo "Usage: $SCRIPT snap_output_dir main_dir"
	exit 1
fi

snap_files=$1
SCRIPTS=$2
. $SCRIPTS/configs.cf 

cd $snap_files
cat $REF_DIR/$REF_NAME.size | while read line
do
	chr=`echo $line | awk '{print $1}'`
	less $REF_GFF | grep -w $chr | grep -e gene -e CDS > $snap_files/$chr.gff	
	$SNAP/gff2zff.pl -sp=$chr $snap_files/$chr.gff
	cat $snap_files/$chr.ann >> $snap_files/$REF_NAME.ann #recommend to do this manually
#	rm $snap_files/$chr.gff $snap_files/$chr.ann
done

$SNAP/fathom $snap_files/$REF_NAME.ann $snap_files/$REF_NAME.dna -gene-stats
$SNAP/fathom $snap_files/$REF_NAME.ann $snap_files/$REF_NAME.dna -validate
$SNAP/fathom $snap_files/$REF_NAME.ann $snap_files/$REF_NAME.dna -categorize 1000
$SNAP/fathom $snap_files/uni.ann $snap_files/uni.dna -export 1000 -plus

mkdir -p params
cd params
$SNAP/forge ../export.ann ../export.dna
