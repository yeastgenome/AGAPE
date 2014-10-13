#ifndef UTIL_EXONS_H
#define UTIL_EXONS_H
#include "main.h"

#define ALL_MODE 0
#define CHR_MODE 1

#define POS_BASE 2
#define LEN_BASE 3

#define SHORT_LEN_TH 101
#define NUM_CHARS 1000
#define NUM_LOOPS 10

void init_exons(struct exons_list *exons, int from, int to);
void assign_gff_exons(FILE *f, struct exons_list *exons, int num_exons);
void assign_gff_exons_chr(FILE *f, struct exons_list *exons, int num_exons, char *chr);
void assign_gff_exons_mode(FILE *f, struct exons_list *exons, int num_exons, char *chr, int mode);
void quick_sort_dec_exons(struct exons_list *a, int lo, int hi, int mode);
void quick_sort_inc_exons(struct exons_list *a, int lo, int hi, int mode);
struct exons_list assign_exons(struct exons_list a);
int quick_search_close_exons(struct exons_list *sorted, int i, int j, int query);
void selection_sort_exons(struct exons_list *a, int num);
int remove_redundant_intervals(struct exons_list *a, int num);
void print_exons_list(struct exons_list *a, int num);

#endif /* UTIL_EXONS_H */
