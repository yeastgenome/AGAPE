#!/bin/sh -l

cfg_all_fungi_dir=/srv/gs1/projects/cherry/giltae/strains/cfg_files/all_fungi
utils_dir=/srv/gs1/projects/cherry/giltae/apps/utils.d
snap_dir=/srv/gs1/projects/cherry/giltae/snap_files
protein_db=/srv/gs1/projects/cherry/giltae/strains/cfg_files/maker/yeast_protein.fasta
Fungi_protein_db=/srv/gs1/projects/cherry/giltae/protein_db/Fungi.protein.fasta
#fasta_dir=/srv/gs1/projects/cherry/giltae/strains/abyss_assembly
annot_ref=/srv/gs1/projects/cherry/giltae/strains/annot/reference
annot_maker=/srv/gs1/projects/cherry/giltae/strains/annot/maker
ref_dir=/srv/gs1/projects/cherry/giltae/RACA/sacCer3.seq.d

gff_dir=/srv/gs1/projects/cherry/giltae/strains/annot/all_annot
#final_gff_dir=/srv/gs1/projects/cherry/giltae/strains/annot/final_annot

genemark_mod=/srv/gs1/projects/cherry/giltae/strains/GeneMark_model/GeneMark_hmm.mod
GeneMark=/srv/gs1/projects/cherry/apps/GeneMark_S/genemark_suite_linux_64/gmsuite/gmsn.pl

strain_name=$1
out_dir=$2
fasta_dir=$3

rm -rf $out_dir/*
mkdir -p $out_dir/temp
fasta_file=$out_dir/$strain_name.fasta

ln -s $fasta_dir/$strain_name.scf.fasta $fasta_file

cur_gff=$gff_dir/$strain_name/$strain_name.genes.gff
#mkdir -p $final_gff_dir/$strain_name
#out_dir=$final_gff_dir/$strain_name
#final_gfff=$final_gff_dir/$strain_name/$strain_name.genes.gff
final_gfff=$out_dir/$strain_name/$strain_name.genes.gff
temp_dir=$out_dir/temp

less $fasta_file | grep ">" | sed 's/>//g'  > $out_dir/scf.list
while read scf_line
do
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
  if [ $scf_len -gt 300 ]
	then
		$utils_dir/pull_fasta_scaf $fasta_file $scf_name > $temp_dir/$scf_name.fasta
		head -1 $temp_dir/$scf_name.fasta > $temp_dir/head.txt
		less $cur_gff | grep -w 'gene' | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/$scf_name.known.genes.gff
		$utils_dir/non_ref $temp_dir/head.txt $temp_dir/$scf_name.known.genes.gff > $temp_dir/$scf_name.non_orf.codex

		num_len=`less "$temp_dir"/"$scf_name".non_orf.codex | wc -l | awk '{print $1}'`

		if [ $num_len -gt 0 ]
		then
			$utils_dir/pull_c $temp_dir/$scf_name.fasta $temp_dir/$scf_name.non_orf.codex >> $out_dir/non_orf.fasta
		fi
	fi
done < $out_dir/scf.list

cd $out_dir
$GeneMark --format=GFF --imod $genemark_mod $out_dir/non_orf.fasta
less $out_dir/non_orf.fasta.gff | sed '/^#/ d' | sed '/^$/d' | awk '{print $1" maker gene "$5" "$6" . "$8" . UNDEF"}' > $out_dir/additional_orfs.gff

mkdir -p $out_dir/blast_out
mkdir -p $out_dir/gff

cd $out_dir

ln -s $Fungi_protein_db $out_dir/Fungi_protein.fasta
makeblastdb -in $out_dir/Fungi_protein.fasta -dbtype prot -parse_seqids -out fungi

echo "#" > $out_dir/temp.codex

blast_out=$out_dir/blast_out
gff=$out_dir/gff

less $out_dir/non_orf.fasta | grep ">" > $out_dir/head.txt
$utils_dir/conv_scf_pos $out_dir/head.txt $out_dir/additional_orfs.gff > $out_dir/temp.genes.gff

less $fasta_file | grep ">" | sed 's/>//g'  > $out_dir/scf.list
while read scf_line 
do
	scf_name=`echo $scf_line | awk '{print $1}'`
	scf_len=`echo $scf_line | awk '{print $2}'`
	
	if [ $scf_len -gt 300 ]
	then
		cd $temp_dir
		less $cur_gff | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/$scf_name.known.genes.gff
		less $out_dir/temp.genes.gff | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/$scf_name.temp.genes.gff
		$utils_dir/merge_gff $out_dir/temp.codex $temp_dir/$scf_name.temp.genes.gff > $temp_dir/$scf_name.non_ref.genes.gff
		$utils_dir/gff2codex $temp_dir/$scf_name.non_ref.genes.gff CDS > $temp_dir/$scf_name.temp.genes.codex
		$utils_dir/reverse_exon_order $temp_dir/$scf_name.temp.genes.codex > $temp_dir/$scf_name.non_ref.genes.codex
		$utils_dir/pull_fasta_scaf $fasta_file $scf_name > $temp_dir/$scf_name.fasta
		$utils_dir/pull_c $temp_dir/$scf_name.fasta $temp_dir/$scf_name.non_ref.genes.codex > $temp_dir/$scf_name.non_ref.genes.dna
		$utils_dir/dna2aa -v $temp_dir/$scf_name.non_ref.genes.dna 1 > $temp_dir/$scf_name.non_ref.genes.aa
		$utils_dir/check_aa $temp_dir/$scf_name.non_ref.genes.gff $temp_dir/$scf_name.non_ref.genes.aa > $temp_dir/$scf_name.filtered.non_ref.genes.gff
		less $temp_dir/$scf_name.non_ref.genes.gff | grep -w 'gene\|match' | awk -v SCF=${scf_name} '{if ($1 == SCF) print $0}' > $temp_dir/$scf_name.non_ref.gff
		$utils_dir/check_start_stop_codons $temp_dir/$scf_name.fasta $temp_dir/$scf_name.non_ref.gff >> $temp_dir/$scf_name.filtered.non_ref.genes.gff
		$utils_dir/check_splice_signal $temp_dir/$scf_name.fasta $temp_dir/$scf_name.filtered.non_ref.genes.gff > $temp_dir/$scf_name.temp.filtered.non_ref.genes.gff
		$utils_dir/merge_gff $out_dir/temp.codex $temp_dir/$scf_name.temp.filtered.non_ref.genes.gff > $temp_dir/all.non_ref.gff
		$utils_dir/gff2codex $temp_dir/all.non_ref.gff > $temp_dir/$scf_name.all.non_ref.codex
		num_len=`less "$temp_dir"/"$scf_name".all.non_ref.codex | wc -l | awk '{print $1}'`
		if [ $num_len -gt 0 ]
		then
			$utils_dir/pull_c $temp_dir/$scf_name.fasta $temp_dir/$scf_name.all.non_ref.codex > $temp_dir/$scf_name.non_ref.dna
			cd $out_dir
		  blastx -db fungi -query $temp_dir/$scf_name.non_ref.dna -outfmt "7 sallacc pident evalue stitle" -out $temp_dir/$scf_name.non_ref.blastx.out
			cd $temp_dir
		  $utils_dir/update_gff_blastx $temp_dir/all.non_ref.gff $temp_dir/$scf_name.non_ref.blastx.out ENSEMBL 80 > $temp_dir/$scf_name.all.non_ref.genes.gff
		fi
		$utils_dir/rm_redun_gff $temp_dir/$scf_name.all.non_ref.genes.gff > $temp_dir/$scf_name.final.genes.gff
		$utils_dir/filter_gff $temp_dir/$scf_name.final.genes.gff MAKER > $temp_dir/$scf_name.additional.genes.gff
		less $temp_dir/$scf_name.additional.genes.gff >> $temp_dir/$scf_name.known.genes.gff
		$utils_dir/merge_gff $out_dir/temp.codex $temp_dir/$scf_name.known.genes.gff > $gff/$scf_name.genes.gff
		mv $temp_dir/$scf_name.non_ref.genes.gff $gff/$scf_name.non_ref.genes.gff
		mv $temp_dir/repeats.gff $gff/$scf_name.non_ref.repeats.gff
		mv $temp_dir/seq.gff $gff/$scf_name.non_ref.seq.gff
		mv $temp_dir/all.non_ref.gff $gff/$scf_name.all.non_ref.genes.gff
		mv $temp_dir/$scf_name.non_ref.blastx.out $blast_out/$scf_name.non_ref.blastx.out
		less $temp_dir/$scf_name.additional.genes.gff >> $out_dir/$strain_name.genemark.genes.gff
		less $gff/$scf_name.genes.gff >> $out_dir/$strain_name.genes.gff
	fi
	cd $out_dir
done < $out_dir/scf.list
rm -rf $temp_dir/*
