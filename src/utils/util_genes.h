#ifndef UTIL_GENES_H
#define UTIL_GENES_H

#define DOUBLE_ORFS 2

#define BLAST_PID_TH 80
#define BLAST_LEN_TH 70

#define POS_BASE 1
#define PID_BASE 2
#define LEN_BASE 3

#define STATS 0
#define SAME_STOP 1
#define OVERLAP 2
#define MISSED 3
#define NOVEL 4
#define IDENTICAL 5
#define CORRECT 6
#define COMPLETE_MISS 7

int get_len_diff(char *gname, struct I reg, bool is_first);
int get_column_in_gff_nth_column(char *column, int nth, char *res1, char *res2);
double get_evalue(char *res1, char *res2, int num_values, bool *is_first); 
float get_pid(char *res1, char *res2, int num_values);
void initialize_scaffolds(struct scaffold *a, int num);
void initialize_genes(struct g_list *a, int num);
void initialize_exons(struct exons_list *a, int num);
void initialize_blast_list(struct blast_list *a, int num);
void quick_sort_dec_genes(struct g_list *a, int lo, int hi, int mode);
void quick_sort_inc_genes(struct g_list *a, int lo, int hi, int mode);
int quick_search_close_genes(struct g_list *sorted, int i, int j, int query);
void assign_genes(struct g_list *new_exons, int loc, struct g_list a);
int input_genes(FILE *f, struct g_list *genes, struct exons_list *exons);
int input_genes_in_gff(FILE *f, struct g_list *genes, struct exons_list *exons);
int input_genes_blastx(FILE *f, struct blast_list *genes, int type, int cutoff_value, struct g_list *gff_genes, int num_gff_genes, struct exons_list *exons);
int count_genes_blastx(FILE *f);
int count_genes(FILE *f, int *num_exons);
int count_genes_in_gff(FILE *f, int *num_exons);
int rm_redun_genes(struct g_list *genes, int from, int to);
void replace_gene_in_list(struct g_list *genes, int cur_id, int tmp_id);
int count_exons(struct g_list *genes, int from, int to);
void assign_exons(struct exons_list *new_exons, int loc, struct exons_list a);
void assign_blast_list(struct blast_list *new_a, int loc, struct blast_list a);
void write_in_gff(struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons);
int rm_overlap_genes(struct g_list *genes, struct exons_list *exons, int from, int to);
void write_orfs_in_other_splicing(struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons, int cut_off);
void print_item_gff(struct g_list gene, struct exons_list *exons);

#endif /* UTIL_GENES_H */
