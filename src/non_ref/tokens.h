#ifndef TOKENS_H
#define TOKENS_H

void parse_last_column_gff(char *line, char *gname, char *seq_name, char *chr, int *b, int *e, float *pid);
void parse_blastx_line(char *line, char *info, int *b, int *e, int type);
int concat_tokens_on_bracket(char *line, int loc, char *index, char *name);
void parse_line(char *line);
void get_start_and_end(char *line, int *b, int *e);
void print_gene(char str, char *gname, char *line);
void print_exons(char *line);
void split_b_and_e(char *line, int *b, int *e);
void rm_temp_num(char *line, char *name);
double AtoF(char *s);

#endif /* TOKENS_H */
