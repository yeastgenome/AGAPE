#!/bin/sh -l

seq_name=$1
annot_dir=$2
ref_dir=$3
ref_name=$4
chain_main=$5
BIN=$6

seq_file=$annot_dir/$seq_name.fasta
chain_dir=$chain_main/$ref_name/$seq_name/chain
temp_dir=$annot_dir/homologs.d/temp

mkdir -p $annot_dir/homologs.d
mkdir -p $annot_dir/homologs.d/intervals
mkdir -p $annot_dir/homologs.d/maf
mkdir -p $annot_dir/homologs.d/fasta
mkdir -p $temp_dir

mkdir -p $ref_dir/intervals.d
ln -s $ref_dir/$ref_name.size $annot_dir/$ref_name.size

while read line
do
	echo $line
	chr_name=`echo $line | awk '{print $1}'`
	e=`echo $line | awk '{print $2}'`
	echo "$chr_name 1 $e" > $ref_dir/intervals.d/$chr_name.interval

	if [ ! -e $chain_dir/$chr_name.chain ]
	then
		echo "$chain_dir/$chr_name.chain not exist or moved to other place"	
	fi
	$BIN/homologs $chain_dir/$chr_name.chain $ref_dir/intervals.d/$chr_name.interval > $annot_dir/homologs.d/intervals/$chr_name.homologs.intervals

	while read line
	do
		first_char=`echo $line | cut -c 1 | awk '{print $1}'`
		if [ "$first_char" != '#' ]
		then
			scaf_name=`echo $line | awk '{print $1}'`
			echo "> $scaf_name" > $temp_dir/cur_scaf
			$BIN/pull_fasta_scaf $seq_file $scaf_name | tail -n +2 >> $temp_dir/cur_scaf
			b=`echo $line | awk '{print $2}'`
			e=`echo $line | awk '{print $3}'`
#			echo "$b, $e"
			b=`expr $b + 1`
			$BIN/lastz T=2 Y=3400 $ref_dir/chr_seq/$chr_name.fa $temp_dir/cur_scaf[$b..$e] --ambiguous=iupac --format=maf >> $annot_dir/homologs.d/maf/$chr_name.homologs.maf
			count=0
			if [ -f $annot_dir/homologs.d/fasta/$chr_name.homologs.fasta ] 
			then
				count=`less "$annot_dir"/homologs.d/fasta/"$chr_name".homologs.fasta | grep -w "$scaf_name" | wc -l`
			fi

			if [ $count -eq 0 ]
			then	
				echo "> $scaf_name" >> $annot_dir/homologs.d/fasta/$chr_name.homologs.fasta
				$BIN/dna $temp_dir/cur_scaf | tail -n +2 >> $annot_dir/homologs.d/fasta/$chr_name.homologs.fasta	
			fi

		fi
	done < $annot_dir/homologs.d/intervals/$chr_name.homologs.intervals
done < $annot_dir/$ref_name.size 
