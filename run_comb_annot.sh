#!/bin/sh

comb_annot=$1
BLAST=$2
out_name=$3
maker_dir=$4
out_dir=$5
snap_dir=$6
SCRIPTS=$7

. $SCRIPTS/configs.cf

cd $comb_annot
rm -rf $comb_annot/ref.*
ln -s $PROTEIN1 $comb_annot/ref1_protein.fasta
$BLAST/makeblastdb -in ref1_protein.fasta -dbtype prot -parse_seqids -out ref
$SCRIPTS/combined_annot.sh $out_name $comb_annot $maker_dir/seq.fasta $out_dir/annot $maker_dir $snap_dir $SCRIPTS SGD 90 $comb_annot # results in $comb_annot/gff/$out_name.genes.gff

rm -rf $comb_annot/ref.*
ln -s $PROTEIN2 $comb_annot/ref_protein.fasta
ln -s $EST2 $comb_annot/ref_est.fasta
ln -s $REPEAT_PROTEIN $comb_annot/te_protein.fasta
ln -s $CFG_DIR/maker_opts.ctl $comb_annot/maker_opts.ctl
ln -s $CFG_DIR/maker_bopts.ctl $comb_annot/maker_bopts.ctl
ln -s $CFG_DIR/maker_exe.ctl $comb_annot/maker_exe.ctl

$BLAST/makeblastdb -in ref_protein.fasta -dbtype prot -parse_seqids -out ref
$SCRIPTS/unannot_regions.sh $maker_dir/seq.fasta $comb_annot/gff/$out_name.genes.gff $comb_annot $SCRIPTS # resutls in $comb_annot/non_orf.fasta
rm $comb_annot/seq.fasta
ln -s $comb_annot/non_orf.fasta $comb_annot/seq.fasta
$SCRIPTS/maker.sh $comb_annot $snap_dir $SCRIPTS # results in $comb_annot/genes.gff
less $comb_annot/non_orf.fasta | grep ">" > $comb_annot/head.txt
mv $comb_annot/genes.gff $comb_annot/add1.genes.gff
$BIN/conv_scf_pos $comb_annot/head.txt $comb_annot/add1.genes.gff > $comb_annot/genes.gff

echo "#" > $comb_annot/$out_name.codex
mkdir -p $comb_annot/more_annot
$SCRIPTS/combined_annot.sh $out_name $comb_annot/more_annot $maker_dir/seq.fasta $comb_annot $comb_annot $snap_dir $SCRIPTS ENSEMBL 80 $comb_annot # results in $comb_annot/more_annot/gff/$out_name.genes.gff
less $comb_annot/more_annot/gff/$out_name.genes.gff >> $comb_annot/gff/$out_name.genes.gff
less $comb_annot/more_annot/blast_out/$out_name.blastx.out >> $comb_annot/blast_out/$out_name.blastx.out
rm -rf $comb_annot/more_annot

rm -rf $comb_annot/seq.fasta
rm -rf $comb_annot/non_orf.fasta
rm -rf $comb_annot/genes.gff
$SCRIPTS/unannot_regions.sh $maker_dir/seq.fasta $comb_annot/gff/$out_name.genes.gff $comb_annot $SCRIPTS # resutls in $comb_annot/non_orf.fasta
$GeneMark  --format=GFF --imod $genemark_mod $comb_annot/non_orf.fasta
less $comb_annot/non_orf.fasta.gff | sed '/^#/ d' | sed '/^$/d' | awk '{print $1" maker gene "$5" "$6" . "$8" . UNDEF"}' > $comb_annot/genes.gff
less $comb_annot/non_orf.fasta | grep ">" > $comb_annot/head.txt
mv $comb_annot/genes.gff $comb_annot/add2.genes.gff
$BIN/conv_scf_pos $comb_annot/head.txt $comb_annot/add2.genes.gff > $comb_annot/genes.gff

echo "#" > $comb_annot/$out_name.codex
mkdir -p $comb_annot/more_annot
$SCRIPTS/combined_annot.sh $out_name $comb_annot/more_annot $maker_dir/seq.fasta $comb_annot $comb_annot $snap_dir $SCRIPTS ENSEMBL 80 $comb_annot # results in $comb_annot/more_annot/gff/$out_name.genes.gff
less $comb_annot/more_annot/gff/$out_name.genes.gff >> $comb_annot/gff/$out_name.genes.gff
less $comb_annot/more_annot/blast_out/$out_name.blastx.out >> $comb_annot/blast_out/$out_name.blastx.out
rm -rf $comb_annot/more_annot
