#!/bin/sh -l

ref_dir=/srv/gs1/projects/cherry/sacCer3
name=$1
seq_dir=$2
out_dir=$3

if [ $# -ne 3 ] && [ $# -ne 4 ]
then
  echo "Usage:  $SCRIPT strain_name seq_dir output_dir [32 or 64]"
  exit 1
elif [ $# == 4 ]
then
	qual_encode=$4
else	
  qual_encode=32
fi

if [ $qual_encode -ne 32 ] && [ $qual_encode -ne 64 ] 
then
  echo "Usage:  $SCRIPT strain_name seq_dir output_dir [32 or 64]"
  exit 1
fi

prep.sh $seq_dir/"$name"_1.fastq $out_dir/"$name"_1.fastq
prep.sh $seq_dir/"$name"_2.fastq $out_dir/"$name"_2.fastq

if [ $qual_encode == 64 ]
then
	bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 -I $ref_dir/sacCer3.fa "$seq_dir/$name"_1.fastq > "$out_dir/$name"_1.fastq.sai
 bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 -I $ref_dir/sacCer3.fa "$seq_dir/$name"_2.fastq > "$out_dir/$name"_2.fastq.sai
else
	bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 $ref_dir/sacCer3.fa "$seq_dir/$name"_1.fastq > "$out_dir/$name"_1.fastq.sai
	bwa aln -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6 -t 1 $ref_dir/sacCer3.fa "$seq_dir/$name"_2.fastq > "$out_dir/$name"_2.fastq.sai
fi

aln_sam.sh $out_dir/$name.bam "$out_dir/$name"_1.fastq "$out_dir/$name"_2.fastq "@RG\tID:ILLUMINA_"$name"\tLB:3\tSM:0167\tPL:illumina"

sam_sort.sh $out_dir/$name.bam $out_dir/$name.sorted.bam 12
sam_index.sh $out_dir/$name.sorted.bam
sam_rm.sh $out_dir/$name.bam

while read line
do 
	chr=`echo $line | awk '{print $1}'`
	bin_sam.sh $chr $out_dir/$chr.bam $out_dir/$name.sorted.bam
	sam_index.sh $out_dir/$chr.bam
	clean_nodup.sh $out_dir/$chr.bam $out_dir/$chr.nodup.bam
	sam_index.sh $out_dir/$chr.nodup.bam 
	rm -rf $out_dir/$chr.bam $out_dir/$chr.bam.bai
	clean_realn.sh $out_dir/$chr.nodup.bam $out_dir/$chr.realn.bam
	sam_index.sh $out_dir/$chr.realn.bam
	rm -rf $out_dir/$chr.nodup.bam $out_dir/$chr.nodup.bam.bai
	clean_recal.sh $out_dir/$chr.realn.bam $out_dir/$chr.recal.temp.bam
	java -Xms5g -Xmx5g -jar $PICARD/AddOrReplaceReadGroups.jar INPUT=$out_dir/$chr.recal.temp.bam OUTPUT=$out_dir/$chr.recal.bam RGID=ILLUMINA_$name RGPL=illumina RGLB=1 RGPU=1111 RGSM=$name VALIDATION_STRINGENCY=LENIENT
	sam_index.sh $out_dir/$chr.recal.bam
	rm -rf $out_dir/$chr.realn.bam $out_dir/$chr.realn.bam.bai $out_dir/$chr.realn.intervals $out_dir/$chr.recal.temp.bam
done < $ref_dir/sacCer3.fa.fai

snp_vcfs=''
indel_vcfs=''
pileup_vcfs=''
cat $ref_dir/sacCer3.fa.fai | ( while read line
do 
	chr=`echo $line | awk '{print $1}'`
	var_snp.sh $out_dir/$chr.snp.gatk.vcf $out_dir/$chr.recal.bam
	var_filter.sh $out_dir/$chr.snp.gatk.vcf
	var_indel.sh $out_dir/$chr.indel.gatk.vcf $out_dir/$chr.recal.bam
	var_filter.sh $out_dir/$chr.indel.gatk.vcf
	var_pileup.sh $out_dir/$chr.pileup.vcf $out_dir/$chr.recal.bam
	snp_vcfs="$snp_vcfs $out_dir/$chr.snp.gatk.vcf"
	indel_vcfs="$indel_vcfs $out_dir/$chr.indel.gatk.vcf"
	pileup_vcfs="$pileup_vcfs $out_dir/$chr.pileup.gatk.vcf"
done

concat_vcf.sh $out_dir/$name.snp.gatk.raw.vcf $snp_vcfs
concat_vcf.sh $out_dir/$name.indel.gatk.raw.vcf $indel_vcfs
concat_vcf.sh $out_dir/$name.pileup.raw.vcf $pileup_vcfs
#merge_vcf.sh $out_dir/$name.snpindel.vcf $out_dir/$name.snp.gatk.raw.vcf $out_dir/$name.indel.gatk.raw.vcf $out_dir/$name.pileup.raw.vcf
annotate.py $out_dir/$name.snp.tsv $out_dir/$name.snp.gatk.raw.vcf
annotate.py $out_dir/$name.indel.tsv $out_dir/$name.indel.gatk.raw.vcf
annotate.py $out_dir/$name.pileup.tsv $out_dir/$name.pileup.raw.vcf

gffs=''
chr=`echo $line | awk '{print $1}'`
var_sv_rpm.sh $out_dir/$name.rpm.gff $out_dir/$name.sorted.bam
var_sv_sra.sh $out_dir/$name.sra.gff $out_dir/$name.sorted.bam
var_sv_rda.sh $out_dir/$name.rda.gff $out_dir/$name.sorted.bam
var_sv_jct.sh $out_dir/$name.jct.gff $out_dir/$name.sorted.bam
var_sv_ctx.sh $out_dir/$name.ctx.gff $out_dir/$name.sorted.bam
gffs="$gffs $out_dir/$name.rpm.gff $out_dir/$name.sra.gff $out_dir/$name.rda.gff $out_dir/$name.jct.gff $out_dir/$name.ctx.gff"

merge_gff.sh $out_dir/$name.svcnv.gff $gffs
annotate.py $out_dir/$name.svcnv.tsv $out_dir/$name.svcnv.raw.gff )
