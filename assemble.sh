#!/bin/sh -l

SCRIPTS=$4
. $SCRIPTS/configs.cf 

cur_dir=$1
seq_name=$2
mode=$3
MIN_CONTIG_LENGTH=100
MIN_PAIRS=5
MINA_NUM=95
K_MER=41

temp_dir=$cur_dir/"$seq_name"_assembly
mkdir -p $temp_dir

if [ ! -e "$ABYSS/abyss-pe" ]
then
	echo "Error: $abyss not exist - maybe not installed or path not edited"
	exit 1
fi

if [ ! -e "$SGA_src/sga-align" ]
then
	echo "Error: $SGA_src/sga-align not exist - maybe not installed or path not edited"
	exit 1
fi

cd $temp_dir
if [ $mode == 1 ]
then
	fname1=$cur_dir/"$seq_name"_1.se.fastq
	lname1=$cur_dir/seq1.fastq
	ln -s $fname1 $lname1
	$ABYSS/abyss-pe aligner=map k=$K_MER name=$seq_name se='seq1.fastq'
elif [ $mode == 2 ]
then
	fname1=$cur_dir/"$seq_name"_1.pe.fastq
	fname2=$cur_dir/"$seq_name"_2.pe.fastq
	fname3=$cur_dir/"$seq_name"_1.se.fastq
	fname4=$cur_dir/"$seq_name"_2.se.fastq
	lname1=$temp_dir/seq1.fastq
	lname2=$temp_dir/seq2.fastq
	lname3=$temp_dir/seq_se1.fastq
	lname4=$temp_dir/seq_se2.fastq
	ln -s $fname1 $lname1
	ln -s $fname2 $lname2
	ln -s $fname3 $lname3
	ln -s $fname4 $lname4
	$ABYSS/abyss-pe aligner=map k=$K_MER name=$seq_name lib='pe1' pe1='seq1.fastq seq2.fastq' se='seq_se1.fastq seq_se2.fastq'
elif [ $mode == 3 ] ## for unmapped reads
then
	fname1=$cur_dir/all_unmapped_reads1.fastq
	lname1=$temp_dir/seq1.fastq
	ln -s $fname1 $lname1
	$ABYSS/abyss-pe aligner=map k=$K_MER name=$seq_name se='seq1.fastq'
elif [ $mode == 4 ] ## for unmapped reads
then
	fname1=$cur_dir/all_unmapped_reads1.fastq
	fname2=$cur_dir/all_unmapped_reads2.fastq
	lname1=$temp_dir/seq1.fastq
	lname2=$temp_dir/seq2.fastq
	ln -s $fname1 $lname1
	ln -s $fname2 $lname2
	$ABYSS/abyss-pe aligner=map k=$K_MER name=$seq_name lib='pe1' pe1='seq1.fastq seq2.fastq'
else
	echo "$mode : unsupported mode (should be 0 or 1)"
	exit 1
fi

chmod 755 $temp_dir/$seq_name-contigs.fa
ln -s $temp_dir/$seq_name-contigs.fa $temp_dir/$seq_name.ctg.fasta

$SGA_src/sga-align --name $seq_name.frag $temp_dir/$seq_name.ctg.fasta $lname1 $lname2
$SGA/sga-bam2de.pl -n $MIN_PAIRS -m $MIN_CONTIG_LENGTH --mina $MINA_NUM --prefix $seq_name.frag $temp_dir/$seq_name.frag.bam
$SGA/sga-astat.py -m $MIN_CONTIG_LENGTH $temp_dir/$seq_name.frag.refsort.bam > $temp_dir/$seq_name.ctg.astat
$SGA/sga scaffold -m $MIN_CONTIG_LENGTH --pe $temp_dir/$seq_name.frag.de -a $temp_dir/$seq_name.ctg.astat -o $temp_dir/$seq_name.scaf \
	$temp_dir/$seq_name.ctg.fasta
$SGA/sga scaffold2fasta -m $MIN_CONTIG_LENGTH -f $temp_dir/$seq_name.ctg.fasta -o $temp_dir/$seq_name.scf.fasta $temp_dir/$seq_name.scaf --write-unplaced --use-overlap

cp $seq_name-contigs.fa $cur_dir/$seq_name.ctg.fasta
cp $seq_name.scf.fasta $cur_dir/

cd $cur_dir
#rm -rf $temp_dir
