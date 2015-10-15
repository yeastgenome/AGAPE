#include "main.h"
#include "util.h"
#include "util_genes.h"

#define NON_GENIC_LEN_TH 100

int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000], scf_name[100];
//	int sum = 0;
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0;
	int num_exons = 0;	
	int i = 0, len = 0;
	int count = 1;

	if( argc != 3 ) {
		fatal("non-ref scaffold-head.txt genes.gff\n");
	}
	
	strcpy(buf, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		while(fgets(buf, 10000, f)) {
			if( sscanf(buf+1, "%s %d %*s", scf_name, &len) != 2 ) {
				fatalf("unexpected fasta header: %s", buf);
			}	
		}
	}

	fclose(f);

	if((f = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[2]);
	}
	else {
		num_genes = count_genes_in_gff(f, &num_exons);
		if( num_genes > 0 ) {
			genes = (struct g_list *) ckalloc(num_genes * sizeof(struct g_list));	
			if( num_exons < num_genes ) num_exons = num_genes;
			exons = (struct exons_list *) ckalloc(num_exons * sizeof(struct exons_list));	
			initialize_genes(genes, num_genes);
			initialize_exons(exons, num_exons);
		}
	}
	fseek(f, 0, SEEK_SET);

	i = input_genes_in_gff(f, genes, exons);	
	if( i != num_genes ) {
		fatalf("gene counter error in %s\n", argv[2]);
	}

	if( num_genes == 0 ) {
		printf(">1 %d %s-1 \n", len, scf_name);
		printf("1 %d\n", len);
	}
	else {
		for( i = 0; i < (num_genes-1); i++ ) {
			if( i == 0 ) {
				if( ( genes[i].txStart - 1 ) > NON_GENIC_LEN_TH ) {
					printf(">1 %d %s-%d\n", genes[i].txStart-1, scf_name, count);
					printf("1 %d\n", genes[i].txStart-1);
					count++;
				}
			}
			else {
				if( ( genes[i+1].txStart - genes[i].txEnd - 1) > NON_GENIC_LEN_TH ) {
					printf(">%d %d %s-%d\n", genes[i].txEnd+1, genes[i+1].txStart-1, scf_name, count);
					printf("%d %d\n", genes[i].txEnd+1, genes[i+1].txStart-1);
					count++;
				}
			}
		}

		if( (len - genes[i].txEnd) > NON_GENIC_LEN_TH ) {
			printf(">%d %d %s-%d\n", genes[i].txEnd+1, len, scf_name, count);
			printf("%d %d\n", genes[i].txEnd+1, len);
			count++;
		}
	}

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	fclose(f);
	return EXIT_SUCCESS;
}
