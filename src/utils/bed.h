#ifndef BED_H
#define BED_H
#include "main.h"

#define FEATURE_LEN 50
#define CHR_LEN 10

struct bed {
  char chr[CHR_LEN];
	struct I reg;
	char feature[FEATURE_LEN];
	float qual;
	char ref[CHR_LEN];
	char alt[CHR_LEN];
	char pass[FEATURE_LEN];
	char info[MAX_LEN];	
};

struct gd_snp {
	int reads1; // for the ref allele
	int reads2; // for the alt allele
	int num_chr; // 0, 1, or 2
	int qual;
};

void init_bed(struct bed *list, int num_rows);

#endif /* BED_H */
