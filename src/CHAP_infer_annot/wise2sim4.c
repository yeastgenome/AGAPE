#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
	FILE *f;
	int from, to, strand, beg, end;
	char buffer[10000], type[100];
	int exons[20][2], exons_num = 0, i;
 
	if ( argc != 5) {
		printf("wise2sim4 wise-file from to strand\n");
		return 1;
	}
	
	f = fopen(argv[1], "r");
	if(f == NULL) {
		printf("Cannot open file %s\n", argv[1]);
		return 1;
	}
	
	from = atoi(argv[2]);
	to = atoi(argv[3]);

	if(argv[4][0] == '+')
		strand = 1;
	else 
		strand = 0;
	
	if(from > to) {
		printf("Improper position specification: %d, %d\n", from, to);
		return 1;
	}

	fgets(buffer, 10000, f);
	while(fgets(buffer, 10000, f)) {
		if(sscanf(buffer, "%s %d %d", type, &beg, &end) != 3)
			continue;
		
		//printf("%s\n", name);
		
		if(strcmp(type, "Gene") == 0) {
			if(strand == 1)
				printf("> %d %d Gene1\n", beg + from - 1, end + from + 2);
			else
				printf("< %d %d Gene1 (complement)\n", to - end - 2, to - beg + 1);
		}
		else {
			if(strand == 1) {
				exons[exons_num][0] = beg + from - 1;
				exons[exons_num][1] = end + from - 1;
			}
			else {
				exons[exons_num][0] = to - end + 1;
				exons[exons_num][1] = to - beg + 1;
			}
			
			exons_num++;
		}
	}
	
	if(strand == 1) {
		for(i=0; i<exons_num - 1; i++)
			printf("%d %d\n", exons[i][0], exons[i][1]);
		printf("%d %d\n", exons[i][0], exons[i][1] + 3);
	}
	else {
		printf("%d %d\n", exons[exons_num - 1][0] - 3, exons[exons_num - 1][1]);
		for(i=exons_num - 2; i>=0; i--)
			printf("%d %d\n", exons[i][0], exons[i][1]);
	}
							
	fclose(f);

	return 0;
}
