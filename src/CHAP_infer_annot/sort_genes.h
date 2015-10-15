#ifndef SORT_GENES_H
#define SORT_GENES_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include "util.h" 

#define LEN_NAME 100

struct I {
	int lower;
	int upper;
};

struct g_list{
	int gid; // a gene id
	char sname[LEN_NAME];
	char gname[LEN_NAME];
  char strand; // if the strand is '+', the direction is non-reverse
	          // otherwise, the strand is reverse
	int txStart, txEnd; // the location of a gene
	int exonCount;
	int exStart, exEnd;
};

struct slist{
	int id;
	int val;
	int val_red; // the number of reduced local alignments
	int sp_state;
	int add_sp_state;
	bool is_x;
};

struct exons_list{
	int fid; // id in the inital list
	struct I reg;	
	int sp_id;
	int val;
};

#endif /* SORT_GENES_H */
