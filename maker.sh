#!/bin/sh -l

cur_dir=$1
snap_dir=$2
SCRIPTS=$3

. $SCRIPTS/configs.cf 

cd $cur_dir 
rm -rf $cur_dir/seq.*.output
rm -rf $cur_dir/maker*
rm -rf $cur_dir/genome.hmm
#echo "genome=seq.fasta" > $cur_dir/maker_opts.ctl
#tail -n +2 $CFG_DIR/maker_opts.ctl >> $cur_dir/maker_opts.ctl
cp $CFG_DIR/maker_opts.ctl $CFG_DIR/maker_bopts.ctl $CFG_DIR/maker_exe.ctl $cur_dir

$SNAP/hmm-assembler.pl $cur_dir/seq.fasta $snap_dir/params > $cur_dir/genome.hmm
$MAKER/maker
$MAKER/gff3_merge -d $cur_dir/seq.maker.output/seq_master_datastore_index.log -o $cur_dir/seq.gff
less $cur_dir/seq.gff | grep repeat > $cur_dir/repeats.gff
less $cur_dir/seq.gff | grep gene > $cur_dir/genes.gff
