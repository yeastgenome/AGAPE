#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "util.h"

#define TRUE 1
#define FALSE 0

int num_col;

char GeneticCode[4][4][4] = { {{'F', 'F', 'L', 'L'}, {'S', 'S', 'S', 'S'}, {'Y', 'Y', '*', '*'}, {'C', 'C', '*', 'W'}},
							 {{'L', 'L', 'L', 'L'}, {'P', 'P', 'P', 'P'}, {'H', 'H', 'Q', 'Q'}, {'R', 'R', 'R', 'R'}},
							 {{'I', 'I', 'I', 'M'}, {'T', 'T', 'T', 'T'}, {'N', 'N', 'K', 'K'}, {'S', 'S', 'R', 'R'}},
							 {{'V', 'V', 'V', 'V'}, {'A', 'A', 'A', 'A'}, {'D', 'D', 'E', 'E'}, {'G', 'G', 'G', 'G'}} };

void dna2protein(char dna[10000], int frame, int vmode) {
	int index, frame1, frame2, frame3;
	
	index = frame - 1;
	
	while(index <= ((int)strlen(dna)) - 3) {
		if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u')
			frame1 = 0;
		else if(dna[index] == 'C' || dna[index] == 'c')
			frame1 = 1;
		else if(dna[index] == 'A' || dna[index] == 'a')
			frame1 = 2;
		else if(dna[index] == 'G' || dna[index] == 'g')
			frame1 = 3;
		else if(dna[index] == 'N' || dna[index] == 'n')
			frame1 = 4;
		else {
			frame1 = 4;
//			fatal("\nWrong DNA format\n");
		}
		
		index++;
			
		if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u')
			frame2 = 0;
		else if(dna[index] == 'C' || dna[index] == 'c')
			frame2 = 1;
		else if(dna[index] == 'A' || dna[index] == 'a')
			frame2 = 2;
		else if(dna[index] == 'G' || dna[index] == 'g')
			frame2 = 3;
		else if(dna[index] == 'N' || dna[index] == 'n')
			frame2 = 4;
		else {
			frame2 = 4;
//			fatal("\nWrong DNA format\n");
		}
		
		index++;
			
		if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u')
			frame3 = 0;
		else if(dna[index] == 'C' || dna[index] == 'c')
			frame3 = 1;
		else if(dna[index] == 'A' || dna[index] == 'a')
			frame3 = 2;
		else if(dna[index] == 'G' || dna[index] == 'g')
			frame3 = 3;
		else if(dna[index] == 'N' || dna[index] == 'n')
			frame3 = 4;
		else {
			frame3 = 4;
//			fatal("\nWrong DNA format\n");
		}
			
		if( (frame1 == 4) || (frame2 == 4) || (frame3 == 4) ) {
			if( vmode == TRUE ) {
				printf("X");
			}
			else {
				fatal("\n N (or n) not allowed\n");
			}
		}
		else printf("%c", GeneticCode[frame1][frame2][frame3]);
		num_col++;

		if( num_col == 50 )
		{
			printf("\n");
			num_col = 0;
		}
			
		index ++;
	}
	
	printf("\n");
}

int main(int argc, char *argv[]) {
	FILE *f;
	int frame = 1;
	char dna[50000], buffer[50000];
	char s[50000];
	char option[10];
	int vmode = FALSE;

	num_col = 0;
  
	if( argc == 4 ) {
		strcpy(option, argv[1]);
		if( (option[0] == '-') && (option[1] == 'v') ) {
			vmode = TRUE;
		}
		else {
			fatal("dna2aa -v dna-file frame\n");
		}
	}
	else if ( argc != 3) {
		fatal("dna2aa dna-file frame\n");
	}
	
	if( vmode == TRUE ) {
		f = fopen(argv[2], "r");
		if(f == NULL) {
			fatalf("Cannot open file %s\n", argv[2]);
		}
	}
	else {
		f = fopen(argv[1], "r");
		if(f == NULL) {
			fatalf("Cannot open file %s\n", argv[1]);
		}
	}
	
	dna[0] = '\0';
	if( vmode == TRUE ) {
		frame = atoi(argv[3]);
	}
	else {
		frame = atoi(argv[2]);
	}

	if(frame > 3 || frame < 1) {
		fatal("Frame should be 1, 2 or 3\n");
	}
	while(fgets(buffer, 10000, f)) {
		if((buffer[0] == '>') || (buffer[0] == '<') || (buffer[0] == '|')) {
			if(strlen(dna) > 0)
				dna2protein(dna, frame, vmode);			
			dna[0] = '\0';
			printf("%s",buffer);
			num_col = 0;
		}
		else {
			if(strlen(buffer) == 1)
				continue;
			sscanf(buffer, "%s", s);
			strcat(dna, s);
		}
	}
	
	if(strlen(dna) > 0)
		dna2protein(dna, frame, vmode);
			
	fclose(f);

	return EXIT_SUCCESS;
}
