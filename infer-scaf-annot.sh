#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to use annotations in a reference species to estimate those for other
# species, via sequence alignments and AUGUSTUS.
#

SCRIPT=`echo $0 | sed -e 's;.*/;;'`    # script name from command line; path removed for msgs

if [ $# -ne 6 ]
then
	echo "Usage:  $SCRIPT out_dir ref_dir ref_name chr_name $BIN SCRIPTS_dir"
	exit 1
fi

SCRIPTS=$6
. $SCRIPTS/configs.cf

out_dir=$1
TEMP=$out_dir/temp.d                            # assorted temporary files
TMP1=$TEMP/per_species                 # used in loops; removed with each pass
TMP2=$TEMP/per_gene
#--------------------------
# $TEMP/m_temp
# $TEMP/temp_loc
# $TEMP/temp_gene
# $TEMP/temp_dna
# $TEMP/temp_aa
#
# $TMP1/$chr_name.maf
# $TMP1/$chr_name.loc
# $TMP1/$chr_name.loc_bound
# $TMP1/temp.exons
#
# $TMP2/temp_maf
# $TMP2/m_genes
# $TMP2/cur_prot
# $TMP2/p_seq_temp
# $TMP2/p_seq
# $TMP2/gene_loc
# $TMP2/loc_file
# $TMP2/final_loc_file
# $TMP2/other_temp_dna
# $TMP2/other_temp_aa
#--------------------------

ref_dir=$2 # chromosome sequence
ref_name=$3
chr_name=$4
BIN=$5
total_num_genes=0

if [ ! -f $ref_dir/chr_seq/$chr_name.fa ]
then
	echo "Sequence matching chromosome name \"$chr_name\" is not found in directory $ref_dir/chr_seq"
	exit 1
fi

if [ ! -f $ref_dir/$ref_name.temp/$chr_name.m_temp ]
then
	echo "Reference annotation file \"$ref_name.temp/$chr_name.m_temp\" is not found."
	exit 1
fi

rm -rf $TEMP; mkdir -p $TEMP
mkdir -p $out_dir/codex

if [ ! -f $out_dir/codex/$chr_name.codex ]
then
	echo "creating $out_dir/codex/$chr_name.codex"
	echo "#### Resetting $TMP1 for $chr_name"
	rm -f $TMP1/*    # kluge for weird bug in "rm -rf" (e.g. with AFS on Macs)
	rm -rf $TMP1; mkdir -p $TMP1
	$BIN/lastz $out_dir/homologs.d/fasta/$chr_name.homologs.fasta[multi] $ref_dir/$ref_name.temp/$chr_name.temp_gene T=2 Y=3400 --ambiguous=iupac --format=maf > $TMP1/$chr_name.maf
	$BIN/extract_gene_cluster_whole $TMP1/$chr_name.maf > $TMP1/$chr_name.loc
	$BIN/gene_boundaries $TMP1/$chr_name.loc 1 > $TMP1/$chr_name.loc_bound

	num=0
	exec 0<$TMP1/$chr_name.loc_bound
	while read line
	do
		#echo $line
		num=`expr $num + 1`
		TMP2=$TEMP/per_gene
		#echo "#### Resetting $TMP2 for $line"
		rm -f $TMP2/*    # kluge for weird bug in "rm -rf" (e.g. with AFS on Macs)
		rm -rf $TMP2; mkdir -p $TMP2
		scf_name=`echo $line | cut -d " " -f1`
		b=`echo $line | cut -d: -f2 | cut -d " " -f2`
		e=`echo $line | cut -d: -f2 | cut -d " " -f3`
		#echo $scf_name $b $e
		$BIN/pull_fasta_scaf $out_dir/homologs.d/fasta/$chr_name.homologs.fasta $scf_name > $TMP1/cur_seq
		$BIN/lastz "$TMP1/cur_seq[$b,$e]" $ref_dir/$ref_name.temp/$chr_name.temp_gene T=2 Y=3400 --ambiguous=iupac --format=maf > $TMP2/temp_maf
		$BIN/find_match $TMP2/temp_maf > $TMP2/m_genes
		#cat $TMP2/m_genes

		is_done=f
		while read aline
		do
			#echo $aline
			cur_name=`echo $aline | awk '{print $1}'`
			direction=`echo $aline | awk '{print $2}'`
			#echo "$cur_name $direction"
			if [ $direction = '+' ] || [ $direction = '-' ]
			then
				$BIN/pull_one_prot $ref_dir/$ref_name.temp/$chr_name.temp_aa $cur_name > $TMP2/cur_prot
				#cat $TMP2/cur_prot
				num_lines=`cat $TMP2/cur_prot | wc -l`
				num_lines=`expr $num_lines - 1`
				len=`tail -$num_lines $TMP2/cur_prot | wc -c`
				len=`expr $len - $num_lines`
				#echo "len = $len"
			fi
			
			num_seq_lines=`cat $TMP1/cur_seq | wc -l`
			num_seq_lines=`expr $num_seq_lines - 1`
			num_nu=`tail -$num_seq_lines $TMP1/cur_seq | wc -c`
			num_nu=`expr $num_nu - $num_seq_lines`
			
			diff=`expr $e - $b`

			if [ $b -lt 31 ] 
			then
				beg=1
			else
				beg=`expr $b - 30`
			fi
 
			end=`expr $e + 30`
			if [ $end -gt $num_nu ] 
			then
				end=$num_nu
			fi

			if [ $direction = '-' ]
			then
				$BIN/dna $beg,$end $TMP1/cur_seq > $TMP2/p_seq_temp
				$BIN/dna -c $TMP2/p_seq_temp > $TMP2/p_seq
			elif [ $direction = '+' ]
			then
				$BIN/dna $beg,$end $TMP1/cur_seq > $TMP2/p_seq
			fi

			num_lines=0
			cur_lines=0

			if [ "$direction" = "+" ] || [ "$direction" = "-" ]
			then
#				echo $diff $len
				$AUGUSTUS/augustus --species=$AUGUSTUS_REF $TMP2/p_seq > $TMP2/gene_loc
#				echo $beg $end $b $e $cur_name $direction $num_nu
				$BIN/gff2sim4 $TMP2/gene_loc $beg $end $b $e $cur_name $direction $num_nu $scf_name > $TMP2/temp_loc
				$BIN/reverse_exon_order $TMP2/temp_loc > $TMP2/final_loc_file
				
				rm -rf $TMP2/gene_loc
				
				tf=f
				num_lines=`cat $TMP2/final_loc_file | wc -l`
				if [ "$num_lines" -ne 0 ]
				then
					$BIN/pull_c $TMP1/cur_seq $TMP2/final_loc_file > $TMP2/other_temp_dna
					$BIN/dna2aa -v $TMP2/other_temp_dna 1 > $TMP2/other_temp_aa
					#echo "dna2aa completed successfully"
					#cat $TMP2/other_temp_aa
					cur_lines=`cat $TMP2/other_temp_aa | wc -l`
					cur_lines=`expr $cur_lines - 1`
					cur_len=`tail -$cur_lines $TMP2/other_temp_aa | wc -c`
					#echo "num of lines: $cur_lines; num of chars: $cur_len"
					cur_len=`expr $cur_len - $cur_lines`
					comp=`expr $len \* 7 / 10`
					if [ "$comp" -gt 25 ]
					then
						comp=25    # the average length of exons - about 170 bp in human
					fi

					if [ "$cur_len" -lt "$comp" ]
					then
						#echo "short aa seq $cur_len < $comp"
						tf=f
					else
						tf=`$BIN/filter_out $TMP2/other_temp_aa $TMP1/cur_seq $len`
					fi
				fi
			fi

			if [ "$cur_lines" -eq 0 ] || [ "$num_lines" -eq 0 ]
			then
				#echo "no lines"
				tf=f
			fi

			#echo "tf = $tf"
			if [ $is_done = 't' ]
			then
				:    # do nothing
			elif [ $tf = 't' ]
			then
				cat $TMP2/final_loc_file >> $TMP1/$scf_name.temp.exons
				total_num_genes=`expr $total_num_genes + 1`
				is_done=t
			elif [ $tf = 'b' ] || [ $tf = 'M' ] || [ $tf = 'P' ]
			then
				$BIN/ext_loc_info $TMP2/final_loc_file "$tf" >> $TMP1/$scf_name.temp.exons
				total_num_genes=`expr $total_num_genes + 1`
				is_done=t
			fi
			#echo "name direction: $cur_name $direction $tf"
		done < $TMP2/m_genes
	done

	num_scf_exon_files=0
	if [ $total_num_genes -gt 0 ]
		then
		for scf_exons in "$TMP1"/*."temp.exons"
		do
			$BIN/sort_genes $scf_exons >> $out_dir/codex/$chr_name.codex
			num_scf_exon_files=`expr $num_scf_exon_files + 1`
		done
	fi

	if [ $num_scf_exon_files == 0 ]
	then
		echo "" > $out_dir/codex/$chr_name.codex
	fi
fi

echo "#### Cleaning up $TEMP"    # rm's '-v' option may help with debugging
#rm -f $TMP2/*    # kluge for weird bug in 'rm -rf' (e.g. with AFS on Macs)
#rm -rf $TMP2
#rm -f $TMP1/*
#rm -rf $TMP1
#rm -f $TEMP/*
#rm -rf $TEMP
