#ifndef UTIL_I_GEN_H
#define UTIL_I_GEN_H
#include "main.h"

void initialize_I_list(struct I *list, int count);
void initialize_scf_I_list(struct scf_I *list, int count);
void initialize_orf_I_list(struct orf_I *list, int count);
int sort_merge_intervals_and_pid(struct I *regs, int num_regs, struct I *new_regs, int *pid);
int sort_merge_intervals(struct I *regs, int num_regs, struct I *new_regs);
int input_scf_I_list(FILE *f, struct scf_I *match_regions, int count);
int find_match_regions(struct I *cur_regions, int num, struct scf_I *match_regions, int num_match_regions, char *scaf_name);
int input_orf_I_list(FILE *f, struct orf_I *match_regions, int count);

#endif /* UTIL_I_GEN_H */
