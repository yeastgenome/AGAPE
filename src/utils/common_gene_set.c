#include "main.h"
#include "util.h"
#include "util_genes.h"

#define LEN_TH 100

int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000], type[100];
//	int sum = 0;
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0;
	int num_exons = 0;	
	int i = 0, j = 0;
	int count = 0;
	struct n_pair *names; 
	int num_list = 0;

	if( argc != 3 ) {
		fatal("common_gene_set list genes.gff\n");
	}
	
	strcpy(buf, "");
	strcpy(type, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		while(fgets(buf, 10000, f)) {
			num_list++;
		}

		names = (struct n_pair *) ckalloc(num_list * sizeof(struct n_pair));
	}
	fseek(f, 0, SEEK_SET);

	i = 0;
	while(fgets(buf, 10000, f)) {
		if( sscanf(buf, "%s %s", names[i].name1, names[i].name2) != 2 ) 
		{
			fatalf("wrong list format: %s", buf);
		}
		i++;
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
	}

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	fclose(f);
	return EXIT_SUCCESS;
}

