#!/bin/sh

out_dir=$1
seq_name=$2
maker_dir=$3
snap_dir=$4
SCRIPTS=$5

. $SCRIPTS/configs.cf

cd $maker_dir 
rm -rf $maker_dir/*.fasta
rm -rf $maker_dir/*.ctl

ln -s $out_dir/"$seq_name".scf.fasta $maker_dir/seq.fasta
ln -s $PROTEIN1 $maker_dir/ref_protein.fasta
ln -s $EST1 $maker_dir/ref_est.fasta
ln -s $REPEAT_PROTEIN $maker_dir/te_protein.fasta
ln -s $CFG_DIR/maker_opts.ctl $maker_dir/maker_opts.ctl
ln -s $CFG_DIR/maker_bopts.ctl $maker_dir/maker_bopts.ctl
ln -s $CFG_DIR/maker_exe.ctl $maker_dir/maker_exe.ctl

$SCRIPTS/maker.sh $maker_dir $snap_dir $SCRIPTS # results in $maker_dir/genes.gff
