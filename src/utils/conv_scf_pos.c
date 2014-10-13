#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"
#include "tokens.h"

int find_scf_id(char *name, struct scaffold *scf, int num);
int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000];
	char word1[100], word2[100];
//	int sum = 0;
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0;
	int num_exons = 0;	
	int i = 0, j = 0, cur_id = 0;
	int num_scf = 0;
	struct scaffold *scf;

	if( argc != 3 ) {
		fatal("conv_scf_pos scaffold-head.txt genes.gff\n");
	}
	
	strcpy(buf, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		while(fgets(buf, 10000, f)) {
			if( sscanf(buf+1, "%s %s %*s", word1, word2) != 2 ) {
				fatalf("unexpected fasta header: %s", buf);
			}	
			else num_scf++;
		}
	}

	if( num_scf > 0 ) {
		scf = (struct scaffold *) ckalloc(sizeof(struct scaffold) * num_scf);
		initialize_scaffolds(scf, num_scf);
	}

	fseek(f, 0, SEEK_SET);
	
	i = 0;
	while(fgets(buf, 10000, f)) {
		if( sscanf(buf+1, "%s %s %*s", word1, word2) != 2 ) {
			fatalf("unexpected fasta header: %s", buf);
		}	
		else {
			strcpy(scf[i].name, word1);
			split_b_and_e(word2, &scf[i].b, &scf[i].e);
			i++;
		}
		if( i > num_scf ) {
			fatalf("scaffold counting error %d\n", num_scf);
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
	}
	else {
		for( i = 0; i < num_genes; i++ ) {
			cur_id = find_scf_id(genes[i].sname, scf, num_scf);
			genes[i].txStart = genes[i].txStart + scf[cur_id].b - 1;
			genes[i].txEnd = genes[i].txEnd + scf[cur_id].b - 1;
			rm_temp_num(genes[i].sname, word1);
			strcpy(genes[i].sname, word1);
			for( j = genes[i].cdsStart; j <= genes[i].cdsEnd; j++ ) {
				exons[j].reg = assign_I(exons[j].reg.lower + scf[cur_id].b - 1, exons[j].reg.upper + scf[cur_id].b - 1);
			}
		}
	}

	write_in_gff(genes, num_genes, exons, num_exons);
	if( num_scf > 0 ) {
		free(scf);
	}

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	fclose(f);
	return EXIT_SUCCESS;
}

int find_scf_id(char *name, struct scaffold *scf, int num)
{
	int res = -1, i = 0;
	
	for( i = 0; i < num; i++ ) {
		if( strcmp(name, scf[i].name) == 0 ) {
			res = i;
		}
	}

	if( res == -1 ) fatalf("%s not in scaffold list\n", name);
	return(res);
}
