#ifndef READ_MAF_H
#define READ_MAF_H

int cal_pid(char *s, char *t, int ncol);
float cal_pid_maf(char *s, char *t, int ncol);
float cal_pid_maf_beg(char *s, char *t, int beg, int ncol);
void read_maf(char *fname, int mode, struct DotList *algns, int *num_algns, int *size1, int *size2);

#endif /* READ_MAF_H */
