#!/bin/sh

if [ $# -ne 4 ] && [ $# -ne 5 ]
then
	echo "Usage: $SCRIPT out_dir output_name AGAPE_main_path sequence1 [or sequence2 for paired end]"	
elif [ $# -eq 5 ]
then
	seq2=$5
fi
out_dir=$1
strain_name=$2

SCRIPTS=$3
. $SCRIPTS/configs.cf

seq1=$4

cd $out_dir
rm -rf $out_dir/*.fastq.sai

ln -s $seq1 $out_dir/"$strain_name"_1.fastq
$BWA/bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 $REF_FASTA $out_dir/"$strain_name"_1.fastq > $out_dir/"$strain_name"_1.fastq.sai

rm -rf $out_dir/aln.bam
rm -rf $out_dir/aln.sorted.bam

if [ $# -eq 5 ]
then
	ln -s $seq2 $out_dir/"$strain_name"_2.fastq
	$BWA/bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 $REF_FASTA $out_dir/"$strain_name"_2.fastq > $out_dir/"$strain_name"_2.fastq.sai
	$BWA/bwa sampe $REF_FASTA $out_dir/"$strain_name"_1.fastq.sai $out_dir/"$strain_name"_2.fastq.sai $out_dir/"$strain_name"_1.fastq $out_dir/"$strain_name"_2.fastq | $SAMTOOLS/samtools view -bo $out_dir/aln.bam -S -
#else
	$BWA/bwa samse $REF_FASTA $out_dir/"$strain_name"_1.fastq.sai $out_dir/"$strain_name"_1.fastq | $SAMTOOLS/samtools view -bo $out_dir/aln.bam -S -
#fi

$SAMTOOLS/samtools sort $out_dir/aln.bam $out_dir/aln.sorted

if [ $# -eq 5 ]
then
	$SAMTOOLS/samtools view -u -f 4 -F 264 $out_dir/aln.sorted.bam > $out_dir/unmapped_temp1.bam
	$SAMTOOLS/samtools view -u -f 8 -F 260 $out_dir/aln.sorted.bam > $out_dir/unmapped_temp2.bam
	$SAMTOOLS/samtools view -u -f 12 -F 256 $out_dir/aln.sorted.bam > $out_dir/unmapped_temp3.bam
	$SAMTOOLS/samtools merge -u $out_dir/all_unmapped.bam $out_dir/unmapped_temp1.bam $out_dir/unmapped_temp2.bam $out_dir/unmapped_temp3.bam
	$SAMTOOLS/samtools sort $out_dir/all_unmapped.bam $out_dir/all_unmapped.sorted

	$BEDTOOLS/bamToFastq -i $out_dir/all_unmapped.sorted.bam -fq $out_dir/all_unmapped_reads1.fastq -fq2 $out_dir/all_unmapped_reads2.fastq
#	BIN/common_reads $out_dir/all_unmapped_reads1.fastq $seq1 > $out_dir/unmapped_reads1.fastq
#	BIN/common_reads $out_dir/all_unmapped_reads2.fastq $seq2 > $out_dir/unmapped_reads2.fastq
else
	$SAMTOOLS/samtools view -u -f 4 $out_dir/aln.sorted.bam > $out_dir/all_unmapped.bam
	$SAMTOOLS/samtools sort $out_dir/all_unmapped.bam $out_dir/all_unmapped.sorted
	$BEDTOOLS/bamToFastq -i $out_dir/all_unmapped.sorted.bam -fq $out_dir/unmapped_reads1.fastq
fi
