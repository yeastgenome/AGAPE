// bed format: start with 0, but this program converts to start with 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[]) {
	FILE *f;
	char gene_chrom[255];
	int from = 0, to = 1, i = 0;
	int bin = 0, txStart = 0, txEnd = 1, cdsStart = 0, cdsEnd = 1, exonCount = 0, exonStarts[1000], exonEnds[1000], id = 0;
	char name[255], chrom[255], strand = '+', name2[255], cdsStartStat[255], cdsEndStat[255], exonFrames[1000], buffer[10000];
	int is_pseudo = FALSE; 
  
	for( i = 0; i < 1000; i++ ) {
		exonStarts[i] = 0;
		exonEnds[i] = 1;
	}
	
	strcpy(name, "");
	strcpy(chrom, "");
	strcpy(name2, "");
	strcpy(cdsStartStat, "");
	strcpy(cdsEndStat, "");
	strcpy(exonFrames, "");
	strcpy(buffer, "");

	if( argc == 6 ) {
		if( strcmp(argv[5], "pseudo") == 0 ) {
			is_pseudo = TRUE;	
		}
		else {
			fatal("extract_gene ref-file chrom from to (or \"pseudo\")\n");
		}
	}
	else if ( argc != 5 ) {
		fatal("extract_gene ref-file chrom from to\n");
	}
	
	strand = '|';

	f = fopen(argv[1], "r");
	if(f == NULL) {
		fatalf("Cannot open file %s\n", argv[1]);
	}
	
	strcpy(gene_chrom, argv[2]);	
	
	from = atoi(argv[3]);
	to = atoi(argv[4]);
	if(from > to) {
		fatalf("Improper position specification: %d, %d\n", from, to);
	}

	while(fscanf(f, "%d %s %s %c %d %d %d %d %d", &bin, name, chrom, &strand, &txStart, &txEnd, &cdsStart, &cdsEnd, &exonCount) != EOF) {
		for(i=0; i<exonCount; i++) 
			fscanf(f, "%d,", &exonStarts[i]);
		for(i=0; i<exonCount; i++) 
			fscanf(f, "%d,", &exonEnds[i]);
		fgets(buffer, 10000, f);
		if(sscanf(buffer, "%d %s %s %s %s", &id, name2, cdsStartStat, cdsEndStat, exonFrames) != 5)
			continue;
		
		//printf("%s\n", name);
		
		if(strcmp(gene_chrom, chrom) != 0)
			continue;
		if((txEnd < from) || (txStart > to))
			continue;
//		if(strncasecmp(gene_name, name2, strlen(gene_name)) != 0)
//			continue;		

		if(strand == '-')
		{
			if( (cdsEnd - cdsStart) > 1 )
			{
				printf("< %d %d %s (complement)\n", cdsStart + 2, cdsEnd + 1, name2);
			}
			else if( is_pseudo == FALSE ) {
				exonCount = 0;
			}
			else {
				printf("< %d %d %s (complement)\n", txStart + 2, txEnd + 1, name2);
				exonCount = 1;
			}
/*
			else 
			{
				printf("< %d %d %s (complement)\n", txStart + 1, txEnd, name2);
			}
*/
		}
		else if(strand == '+')
		{
			if( (cdsEnd - cdsStart) > 1 )
			{
				printf("> %d %d %s\n", cdsStart + 2, cdsEnd + 1, name2);
			}
			else if( is_pseudo == FALSE ) {
				exonCount = 0;
			}
			else
			{
				printf("> %d %d %s\n", txStart + 2, txEnd + 1, name2);
				exonCount = 1;
			}
		}
/*		if(exonCount > 1) {
			printf("%d %d\n", cdsStart, exonEnds[0]);
			for(i=1; i<exonCount - 1; i++)
				printf("%d %d\n", exonStarts[i], exonEnds[i]);
			printf("%d %d\n", exonStarts[i], cdsEnd);
		}
		else {
			printf("%d %d\n", cdsStart, cdsEnd);
		}*/
		
		if( ( is_pseudo == TRUE ) && (exonCount == 1) ) {
			printf("%d %d\n", txStart + 1, txEnd);
		}
		else {
			if( strand == '-' ) {
				for(i=0; i<exonCount; i++) {
					if((exonStarts[i] > cdsEnd || exonEnds[i] < cdsStart))
						continue;
					if(exonStarts[i] < cdsStart)
						printf("%d ", cdsStart + 2);
					else
						printf("%d ", exonStarts[i] + 2);
					if(exonEnds[i] > cdsEnd + 1)
						printf("%d\n", cdsEnd + 1);
					else
						printf("%d\n", exonEnds[i] + 1);
				}			
			}
			else if( strand == '+' ) {
				for(i=0; i<exonCount; i++) {
					if((exonStarts[i] > cdsEnd || exonEnds[i] < cdsStart))
						continue;
					if(exonStarts[i] < cdsStart)
						printf("%d ", cdsStart + 2);
					else
						printf("%d ", exonStarts[i] + 2);
					if(exonEnds[i] > cdsEnd + 1)
						printf("%d\n", cdsEnd + 1);
					else
						printf("%d\n", exonEnds[i] + 1);
				}
			}
		}
	}
	
	fclose(f);

	return EXIT_SUCCESS;
}
