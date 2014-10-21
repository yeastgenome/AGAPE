#ifndef INPUT_BLASTN_H
#define INPUT_BLASTN_H
#include "main.h"

#define BLASTN_LEN_TH 0.8
#define EVAL_CUTOFF 1e-06
#define PID_CUTOFF 70

bool is_name_already_in(char *name, struct sp_list *list, int num);
int count_genes_blastn(FILE *f, int *max_hits);
int input_genes_blastn(FILE *f, struct sp_list **genes, int pid_cutoff, bool is_print);

#endif /* INPUT_BLASTN_H */
