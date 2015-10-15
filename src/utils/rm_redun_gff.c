#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"

int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000], type[100];
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0;
	int num_exons = 0;	
	int num = 0;
	int i = 0, j = 0;

	if( argc != 2 ) {
		fatal("rm_redun_gff gff\n");
	}
	
	strcpy(buf, "");
	strcpy(type, "");
	
	if((f = ckopen(argv[1], "r")) == NULL )
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
	num = input_genes_in_gff(f, genes, exons);	
	if( num != num_genes ) {
		fatalf("gene counter error in %s\n", argv[1]);
	}

	if( num_genes > 0 ) {
		quick_sort_inc_genes(genes, 0, num_genes-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes ) {
		j = 0;
    while( ((i+j) < num_genes) && (genes[i].txStart == genes[i+j].txStart )) j++;
    quick_sort_dec_genes(genes, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	if( num_genes > 0 ) {
		num_genes = rm_overlap_genes(genes, exons, 0, num_genes-1);
	}

	write_in_gff(genes, num_genes, exons, num_exons);

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	fclose(f);
	return EXIT_SUCCESS;
}
