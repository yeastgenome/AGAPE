#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	char buffer[10000];
	int *c_beg, *c_end, *g_beg, *g_end;
	char **c_name, **g_name;
	int cluster_num = 0, gene_num = 0, i = 0, shift = 0;
  
	if ( argc != 3) {
		printf("gene_boundaries cluster-file shift\n");
		return EXIT_FAILURE;
	}
	
	f = fopen(argv[1], "r");
	if(f == NULL) {
		printf("Cannot open file %s\n", argv[1]);
		return EXIT_FAILURE;
	}
	
	cluster_num = 0;
	while(fgets(buffer, 10000, f)) cluster_num++; 

	if( cluster_num == 0 ) return EXIT_SUCCESS;
	else {
		c_beg = (int *) ckalloc(sizeof(int) * cluster_num);
		c_end = (int *) ckalloc(sizeof(int) * cluster_num);
		g_beg = (int *) ckalloc(sizeof(int) * cluster_num);
		g_end = (int *) ckalloc(sizeof(int) * cluster_num);
		c_name = (char **) ckalloc(sizeof(char *) * cluster_num);
		g_name = (char **) ckalloc(sizeof(char *) * cluster_num);
		for( i = 0; i < cluster_num; i++ ) {
			c_beg[i] = 0;
			c_end[i] = 0;
			g_beg[i] = 0;
			g_end[i] = 0;
			c_name[i] = ckalloc(sizeof(char) * 1000);
			g_name[i] = ckalloc(sizeof(char) * 1000);
			strcpy(c_name[i], "");
			strcpy(g_name[i], "");
		}
	}

	fseek(f, 0, SEEK_SET);

	i = 0;
	while( fgets(buffer, 10000, f) ) {
		sscanf(buffer, "%s %d - %d", c_name[i], &c_beg[i], &c_end[i]);
		i++;
	}
	
	g_beg[0] = c_beg[0];
	g_end[0] = c_end[0];
	strcpy(g_name[0], c_name[0]);
	gene_num = 1;
	
	for(i=1; i<cluster_num; i++) {
		if((strcmp(g_name[gene_num-1], c_name[i]) == 0 ) && (g_beg[gene_num-1] - c_end[i]) < 100) {
			g_beg[gene_num-1] = c_beg[i];
		}
		else {
			strcpy(g_name[gene_num], c_name[i]);	
			g_beg[gene_num] = c_beg[i];
			g_end[gene_num] = c_end[i];
			gene_num++;
		}
	}
			
	fclose(f);
	
	shift = atoi(argv[2]);
	
	for(i=0; i<gene_num; i++) {
		//printf("gene %d : %d\t%d\n", i+1, shift + g_beg[i] - 2, shift + g_end[i] + 2);
		printf("%s %d : %d\t%d\n", g_name[i], i+1, shift + g_beg[i], shift + g_end[i]);
	}
	
	return EXIT_SUCCESS;
}
