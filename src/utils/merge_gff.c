#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"

int comp_g_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, struct g_list *genes2, int num_genes2, struct exons_list *exons2, int num_exons2, struct g_list *genes3, struct exons_list *exons3, int *num_exons3);

int main(int argc, char *argv[]) {
	FILE *f, *g;
//	int b = 0, e = 0;
	char buf[10000], type[100];
//	int sum = 0;
	struct g_list *genes1, *genes2, *genes3;
	struct exons_list *exons1, *exons2, *exons3;
	int num_genes1 = 0, num_genes2 = 0, num_genes3 = 0;
	int num_exons1 = 0, num_exons2 = 0, num_exons3 = 0;	
	int num = 0;
	int i = 0, j = 0;
	char seq_name[100];

	strcpy(seq_name, "");
	if( argc == 4 ) {
		strcpy(seq_name, argv[3]);
	}
	else if( argc != 3 ) {
		fatal("merge_gff codex gff2 (seq_name)\n");
	}
	
	strcpy(buf, "");
	strcpy(type, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		num_genes1 = count_genes(f, &num_exons1);

		if( num_genes1 > 0 ) {
			genes1 = (struct g_list *) ckalloc(num_genes1 * sizeof(struct g_list));	
			exons1 = (struct exons_list *) ckalloc(num_exons1 * sizeof(struct exons_list));	
			if( num_exons1 < num_genes1 ) num_exons1 = num_genes1;
			initialize_genes(genes1, num_genes1);
			initialize_exons(exons1, num_exons1);
		}
	}
	fseek(f, 0, SEEK_SET);
	
	if((g = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[2]);
	}
	else {
		num_genes2 = count_genes_in_gff(g, &num_exons2);
		if( num_genes2 > 0 ) {
			genes2 = (struct g_list *) ckalloc(num_genes2 * sizeof(struct g_list));	
			if( num_exons2 < num_genes2 ) num_exons2 = num_genes2;
			exons2 = (struct exons_list *) ckalloc(num_exons2 * sizeof(struct exons_list));	
			initialize_genes(genes2, num_genes2);
			initialize_exons(exons2, num_exons2);
		}
	}
	fseek(g, 0, SEEK_SET);

	num = input_genes(f, genes1, exons1);	
	if( num != num_genes1 ) {
		fatalf("gene counter error in %s\n", argv[1]);
	}

	num = input_genes_in_gff(g, genes2, exons2);	
	if( num != num_genes2 ) {
		fatalf("gene counter error in %s\n", argv[2]);
	}

	if( num_genes1 > 0 ) {
		quick_sort_inc_genes(genes1, 0, num_genes1-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes1 ) {
		j = 0;
    while( ((i+j) < num_genes1) && (genes1[i].txStart == genes1[i+j].txStart )) j++;
    quick_sort_dec_genes(genes1, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	if( num_genes2 > 0 ) {
		quick_sort_inc_genes(genes2, 0, num_genes2-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes2 ) {
		j = 0;
    while( ((i+j) < num_genes2) && (genes2[i].txStart == genes2[i+j].txStart )) j++;
    quick_sort_dec_genes(genes2, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	if( num_genes1 > 0 ) {
		num_genes1 = rm_redun_genes(genes1, 0, num_genes1-1);
		num_exons3 = count_exons(genes1, 0, num_genes1-1);
	}

	if( num_genes2 > 0 ) {
		num_genes2 = rm_redun_genes(genes2, 0, num_genes2-1);
		num_exons3 = num_exons3 + count_exons(genes2, 0, num_genes2-1);
	}
	num_genes3 = num_genes1 + num_genes2;

	if( num_genes3 > 0 ) {
		genes3 = (struct g_list *) ckalloc(num_genes3 * sizeof(struct g_list));
		if( num_exons3 < num_genes3 ) num_exons3 = num_genes3;
		exons3 = (struct exons_list *) ckalloc(num_exons3 * sizeof(struct exons_list));
		initialize_genes(genes3, num_genes3);
		initialize_exons(exons3, num_exons3);
	}

	num_genes3 = comp_g_lists(genes1, num_genes1, exons1, num_exons1, genes2, num_genes2, exons2, num_exons2, genes3, exons3, &num_exons3);

	if( num_genes3 > 0 ) {
		quick_sort_inc_genes(genes3, 0, num_genes3-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes3 ) {
		j = 0;
    while( ((i+j) < num_genes3) && (genes3[i].txStart == genes3[i+j].txStart )) j++;
    quick_sort_dec_genes(genes3, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	write_in_gff(genes3, num_genes3, exons3, num_exons3);

	if( num_genes1 > 0 ) {
		free(genes1);
		free(exons1);
	}

	if( num_genes2 > 0 ) {
		free(genes2);
		free(exons2);
	}

	if( num_genes3 > 0 ) {
		free(genes3);
		free(exons3);
	}
	fclose(f);
	fclose(g);
	return EXIT_SUCCESS;
}

// gene names in genes1 are known and not in genes2
int comp_g_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, struct g_list *genes2, int num_genes2, struct exons_list *exons2, int num_exons2, struct g_list *genes3, struct exons_list *exons3, int *num_exons3)
{
	int i = 0, j = 0, k = 0;
	int num_exons = 0;
	struct I cur, tmp;

	cur = assign_I(0, 1);
	tmp = assign_I(0, 1);

	for( i = 0; i < num_genes1; i++ ) 
	{
		cur = assign_I(genes1[i].txStart, genes1[i].txEnd);

		if( num_genes2 > 0 ) {
			j = 0;
			tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
		}

		while( (j <= (num_genes2-1)) && (overlap(cur, tmp) == false) ) {
			j++;

			if( j < (num_genes2-1) ) {
				tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			}
		}

		while( (j <= (num_genes2-1)) && (overlap(cur, tmp) == true) ) {
			if( equal(cur, tmp) == true ) {
				genes2[j].type = REDUN;
			}			
			else if( subset(cur, tmp) == true ) {
				genes2[i].type = REDUN;
			}
			j++;
			if( j < (num_genes2-1) ) {
				tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			}
		}
	}

	for( i = 0; i < num_genes1; i++ ) {
		assign_genes(genes3, i, genes1[i]);
		genes3[i].cdsStart = num_exons;
		for( k = genes1[i].cdsStart; k <= genes1[i].cdsEnd; k++ ) {
			if( k >= num_exons1 ) {
				fatalf("%d excceeds the range [0,%d]\n", k, num_exons1);
			}
			assign_exons(exons3, num_exons, exons1[k]);
			num_exons++;
		}
		genes3[i].cdsEnd = num_exons - 1;
		if( genes3[i].cdsEnd < genes3[i].cdsStart ) {
			fatalf("exons assignments error in %d-%d\n", genes3[i].txStart, genes3[i].txEnd);
		}
	}
	
	for( j = 0; j < num_genes2; j++ ) {
		if( genes2[j].type != REDUN ) {
			assign_genes(genes3, i, genes2[j]);
			genes3[i].cdsStart = num_exons;
			for( k = genes2[j].cdsStart; k <= genes2[j].cdsEnd; k++ ) {
				if( k >= num_exons2 ) {
					fatalf("%d excceeds the range [0,%d]\n", k, num_exons2);
				}
				assign_exons(exons3, num_exons, exons2[k]);
				num_exons++;
			}
			genes3[i].cdsEnd = num_exons - 1;
			
			if( genes3[i].cdsEnd < genes3[i].cdsStart ) {
				fatalf("exons assignments error in %d-%d\n", genes3[i].txStart, genes3[i].txEnd);
			}
			i++;
		}
	} 
	*num_exons3 = num_exons;
	return(i);
}
