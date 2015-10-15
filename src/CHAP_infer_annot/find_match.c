#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include "util.h"
#include "util_sort.h"

#define BIG 1000000
#define TRUE 0
#define FALSE 1

int main(int argc, char *argv[]) {
	FILE *fp;
	char S[BIG], T[BIG];
	char *direction, **gene;
	int gene_index = 0, *scores;
	int b = 0, e = 0;
	char *strand;
	char *status;
	int temp_score = 0;
	int num_algns = 0;
	struct slist *st;
	int i = 0, j = 0, k = 0;
	char best_matches[3][100];
	char best_dir[3];
	int is_new = TRUE;
  
	for( i = 0; i < 3; i++ ) {
		strcpy(best_matches[i], "");
		best_dir[i] = '+';
	}

	if ( argc != 2) {
		fatal("find_match maf-format\n");
	}
	
	fp = fopen(argv[1], "r");
	if(fp == NULL) {
		fatalf("Cannot open file %s\n", argv[1]);
	}

	num_algns = 0;
	while(fgets(S, BIG, fp)) {
		if( S[0] == 'a' ) num_algns++;
	}

	strand = (char *) ckalloc(sizeof(char)+1);
	if( num_algns > 0 ) {
		direction = (char *) ckalloc(sizeof(char) * num_algns);
		scores = (int *) ckalloc(sizeof(int) * num_algns);
		gene = (char **) ckalloc(sizeof(char *) * num_algns);
	}
	else {
		direction = (char *) ckalloc(sizeof(char));
		scores = (int *) ckalloc(sizeof(int));
		gene = (char **) ckalloc(sizeof(char *));
	}
	strcpy(strand, "|");
	for( i = 0; i < num_algns; i++ ) {
		gene[i] = (char *) ckalloc(sizeof(char) * 100);		
		direction[i] = '|';
		scores[i] = 0;
		strcpy(gene[i], "");
	}
	
	fseek(fp, 0, SEEK_SET);
	if (fgets(S, BIG, fp) == NULL || strncmp(S, "##maf version", 13)) {
		fatalf("%s is not a maf file\n", S);
	}

	if( num_algns > 0 ) {
		while (S[0] == '#') {
			if(S[0] == '#'){
				while(S[0] == '#') {
					if((status = fgets(S, BIG, fp)) == NULL) {
						fatalf("no alignments in %s\n", argv[1]);
					}
				}
			}
		}	
	}

	gene_index = 0;
	for(i=0; i < num_algns; i++)
		scores[i] = 0;

	while ((num_algns > 0) && (status != NULL) && (strstr(S, "eof") == NULL)) {
		if(strncmp(S, "a ", 2)) {
	  	fatalf("expecting an a-line in %s, saw %s", argv[1], S);
		}
		if(sscanf(S+8, "%d", &temp_score) == 1) scores[gene_index] = temp_score;
		else {
			fatalf("%s is not correct maf format\n", S);
		}

    if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL)) {
      fatalf("cannot find alignment in %s", argv[1]);
		}

		if ((sscanf(T, "%*s %s %d %d %s %*s", gene[gene_index], &b, &e, strand) == 4)) {
			if( (strstr(strand, "-") != NULL) ) direction[gene_index] = '-';
			else if( (strstr(strand, "+") != NULL) ) direction[gene_index] = '+';
		}
		else
		{
			fatalf("%s is not correct maf format\n", T);
		}
		
	  if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n')) {
      fatalf("bad alignment end in %s", argv[1]);
		}
    status = fgets(S, BIG, fp);

		if( direction[gene_index] == '|' ) {}
		else {
			gene_index++;
		}
	}

	fclose(fp);
	
	if( gene_index > 0 ) st = (struct slist *) ckalloc(sizeof(struct slist) * gene_index);

	for( i = 0; i < gene_index; i++ ) {
		st[i].id = i;
		st[i].val = scores[i];
	}
	
	if( gene_index > 0 ) {
		quick_sort_dec(st, 0, gene_index-1);

		i = 0;
		j = 0;
		while( (i < gene_index) && (j < 3) ) {
			for( k = 0; k < j; k++ ) {
				if( strcmp(best_matches[k], gene[st[i].id]) == 0 ) {
					is_new = FALSE;
				}
			}

			if( is_new == TRUE ) {
				strcpy(best_matches[j], gene[st[i].id]);
				best_dir[j] = direction[st[i].id];
				j++;
			}
			i++;
			is_new = TRUE;
		}

		if( j >= 3 ) {
			for( i = 0; i < 3; i++ ) {
				printf("%s %c\n", best_matches[i], best_dir[i]);
			}
		}
		else {
			for( i = 0; i < j; i++ ) {
				printf("%s %c\n", best_matches[i], best_dir[i]);
			}
		}
		free(st);
	}

	free(direction); 
	free(strand);
	free(scores);
	for( i = 0; i < num_algns; i++ ) {
		free(gene[i]);
	}
	free(gene);
	return EXIT_SUCCESS;
}
