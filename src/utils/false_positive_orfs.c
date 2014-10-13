#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"

int main(int argc, char *argv[]) {
	FILE *f;
//	int b = 0, e = 0;
	char buf[10000];
//	int sum = 0;
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0, num_exons = 0;
	int num = 0;
	int cutoff_value = BLAST_PID_TH;

	if( argc == 3 ) {
		cutoff_value = atoi(argv[2]);
	}
	else if( argc != 2 ) {
		fatal("false_positive_orfs gff (cut_off)\n");
	}
	else {
	}
	
	strcpy(buf, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		num_genes = count_genes_in_gff(f, &num_exons);
		if( num_genes > 0 ) {
			genes = (struct g_list *) ckalloc(num_genes * sizeof(struct g_list));	
			exons = (struct exons_list *) ckalloc(num_exons * sizeof(struct exons_list));	
			initialize_genes(genes, num_genes);
			initialize_exons(exons, num_exons);
		}
	}
	fseek(f, 0, SEEK_SET);

	num = input_genes_in_gff(f, genes, exons);	
	if( num != num_genes ) {
		fatalf("gene counter error in %s\n", argv[1]);
	}

	fclose(f);
	
	write_orfs_in_other_splicing(genes, num_genes, exons, num_exons, cutoff_value);

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	return EXIT_SUCCESS;
}
