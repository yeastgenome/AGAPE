#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	int b = 0, e = 0;
	char buf[10000], type[100];
	int sum = 0;
	int num_genes = 0;
 
	if ( argc != 2 ) {
		fatal("sumgff gff_annotation\n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	
	strcpy(buf, "");

	num_genes = 0;
	sum = 0;
	while( fgets(buf, 10000, f)) { 
		if( sscanf(buf, "%*s %*s %s %d %d %*s", type, &b, &e) != 3 ) {
			fatalf("wrong gff format: %s", buf);
		}

		if( (strcmp(type, "gene") == 0) ) {
			sum = sum + (e - b + 1);
			num_genes++;
		}
	}
	
	printf("%d %d\n", num_genes, sum);
	fclose(f);
	return EXIT_SUCCESS;
}
