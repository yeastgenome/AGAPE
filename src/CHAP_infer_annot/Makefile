include ../Flags.mk
DEST = ../../bin

PROGS = get_gff_gene_bound pull_c recal_loc extract_all_gene make_ugname get_gene_bound dna2aa extract_gene_cluster_whole gene_boundaries wise2sim4 find_match ch_gname filter_out pull_one_prot ext_loc_info sort_genes dna gff2temp get_gene_bound

all : $(PROGS)

gff2temp :
	$(CC) $(CFLAGS) util.c util_annot.c gff2temp.c $(LFLAGS) -o $(DEST)/gff2temp

pull_c :
	$(CC) $(CFLAGS) pull_c.c util.c seq.c nib.c charvec.c $(LFLAGS) -o $(DEST)/pull_c

recal_loc :
	$(CC) $(CFLAGS) recal_loc.c $(LFLAGS) -o $(DEST)/recal_loc

extract_all_gene :
	$(CC) $(CFLAGS) extract_all_gene.c util.c $(LFLAGS) -o $(DEST)/extract_all_gene

make_ugname :
	$(CC) $(CFLAGS) make_ugname.c $(LFLAGS) -o $(DEST)/make_ugname

get_gene_bound :
	$(CC) $(CFLAGS) get_gene_bound.c $(LFLAGS) -o $(DEST)/get_gene_bound

get_gff_gene_bound :
	$(CC) $(CFLAGS) get_gff_gene_bound.c $(LFLAGS) -o $(DEST)/get_gff_gene_bound

dna2aa :
	$(CC) $(CFLAGS) dna2aa.c util.c $(LFLAGS) -o $(DEST)/dna2aa

extract_gene_cluster_whole :
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c chain.c mz_scores.c extract_gene_cluster.c $(LFLAGS) -o $(DEST)/extract_gene_cluster_whole

gene_boundaries :
	$(CC) $(CFLAGS) util.c  gene_boundaries.c $(LFLAGS) -o $(DEST)/gene_boundaries

wise2sim4 :
	$(CC) $(CFLAGS) wise2sim4.c $(LFLAGS) -o $(DEST)/wise2sim4

find_match :
	$(CC) $(CFLAGS) find_match.c util.c util_sort.c $(LFLAGS) -o $(DEST)/find_match

ch_gname :
	$(CC) $(CFLAGS) ch_gname.c util.c $(LFLAGS) -o $(DEST)/ch_gname

filter_out :
	$(CC) $(CFLAGS) filter_out.c util.c $(LFLAGS) -o $(DEST)/filter_out

pull_one_prot :
	$(CC) $(CFLAGS) pull_one_prot.c $(LFLAGS) -o $(DEST)/pull_one_prot

ext_loc_info :
	$(CC) $(CFLAGS) ext_loc_info.c $(LFLAGS) -o $(DEST)/ext_loc_info

sort_genes :
	$(CC) $(CFLAGS) sort_genes.c util.c util_i.c util_genes.c $(LFLAGS) -o $(DEST)/sort_genes

dna :
	$(CC) $(CFLAGS) dna.c seq.c util.c nib.c charvec.c $(LFLAGS) -o $(DEST)/dna

