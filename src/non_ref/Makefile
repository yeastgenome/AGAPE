include ./Flags.mk
DEST = ../../bin

PROGS : filter_blast homologous_orfs find_homologous_orfs mask_out_unex_match common_orf_intervals common_scf_intervals merge_scf_intervals scf_lift_over cal_avg_pid 
#find_inserted_intervals 

all : $(PROGS)

filter_blast :
	$(CC) $(CFLAGS) filter_blast.c input_blastn.c tokens.c util.c $(LFLAGS) -o $(DEST)/filter_blast

homologous_orfs :
	$(CC) $(CFLAGS) identify_homologous_orfs_with_blastn.c input_blastn.c tokens.c util.c $(LFLAGS) -o $(DEST)/homologous_orfs

find_homologous_orfs :
	$(CC) $(CFLAGS) identify_homologous_orfs.c read_maf.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/find_homologous_orfs

mask_out_unex_match :
	$(CC) $(CFLAGS) mask_out_unex_match.c util_I_gen.c read_maf.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/mask_out_unex_match
	
common_orf_intervals :
	$(CC) $(CFLAGS) common_orf_intervals.c util_I_gen.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/common_orf_intervals

common_scf_intervals :
	$(CC) $(CFLAGS) common_scf_intervals.c util_I_gen.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/common_scf_intervals

merge_scf_intervals :
	$(CC) $(CFLAGS) merge_scf_intervals.c util_I_gen.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/merge_scf_intervals

scf_lift_over :
	$(CC) $(CFLAGS) scf_lift_over.c read_maf.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/scf_lift_over

find_inserted_intervals :
	$(CC) $(CFLAGS) find_inserted_intervals.c util_I_gen.c read_maf.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/find_inserted_intervals

cal_avg_pid :
	$(CC) $(CFLAGS) cal_avg_pid.c util_I_gen.c read_maf.c kd_tree.c regions.c util_gen.c util_i.c util.c $(LFLAGS) -o $(DEST)/cal_avg_pid
