#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"

void fill_gname(struct g_list *genes, int num_genes, struct blast_list *blastx_res, int num_blastx);

int main(int argc, char *argv[]) {
	FILE *f;
//	int b = 0, e = 0;
	char buf[10000];
	int type = SGD;
//	int sum = 0;
	struct g_list *genes;
	struct exons_list *exons;
	struct blast_list *blastx_res;
	int num_genes = 0, num_exons = 0;
	int num_blastx = 0;
	int num = 0;
	int i = 0, j = 0;
	int cutoff_value = BLAST_PID_TH;

	if( argc == 5 ) {
		cutoff_value = atoi(argv[4]);
	}

	if( argc >= 4 ) {
		if( strcmp(argv[3], "SGD") == 0 ) {
			type = SGD;
		}
		else if( strcmp(argv[3], "ENSEMBL") == 0 ) {
			type = ENSEMBL;
		}
		else {
			type = SGD;
		}
	}
	else if( argc != 3 ) {
		fatal("update_gff_with_blastx_output gff blastx_output\n");
	}
	else {
		type = SGD;
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

	fclose(f);
	
	if((f = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[2]);
	}
	else {
		num_blastx = count_genes_blastx(f);
		if( num_blastx > 0 ) {
			blastx_res = (struct blast_list *) ckalloc(num_blastx * sizeof(struct blast_list));	
			initialize_blast_list(blastx_res, num_blastx);
		}
	}	
	fseek(f, 0, SEEK_SET);
	num_blastx = input_genes_blastx(f, blastx_res, type, cutoff_value, genes, num_genes, exons);	
	
	fill_gname(genes, num_genes, blastx_res, num_blastx);
	write_in_gff(genes, num_genes, exons, num_exons);

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}

	if( num_blastx > 0 ) {
		free(blastx_res);
	}
	fclose(f);
	return EXIT_SUCCESS;
}

void fill_gname(struct g_list *genes, int num_genes, struct blast_list *blastx_res, int num_blastx)
{
	int i = 0, j = 0; 
	struct I cur_reg;
		
	cur_reg = assign_I(0, 1);

	for( i = 0; i < num_genes; i++ ) {
		j = 0;
		cur_reg = assign_I(genes[i].txStart, genes[i].txEnd);
		strcpy(genes[i].gname, "UNDEF");
		strcpy(genes[i].chrname, "");
		strcpy(genes[i].info, "");
		while((j < num_blastx) && (equal(cur_reg, blastx_res[j].reg) == false)) {
			j++;
		} 

		if( j < num_blastx ) {
			if( (int)(blastx_res[j].pid1) != -1 ) {
				strcpy(genes[i].gname, blastx_res[j].name1);
				sprintf(genes[i].chrname, "%s,%s,%.2f",blastx_res[j].name1, blastx_res[j].info1, blastx_res[j].pid1);
			}

			if( (int)(blastx_res[j].pid2) != -1 ) {
				sprintf(genes[i].info, "%s,%s,%.2f",blastx_res[j].name2, blastx_res[j].info2, blastx_res[j].pid2);
				genes[i].ortho_id = DOUBLE_ORFS;
			}
		}
	}
}
