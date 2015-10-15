#include "sort_genes.h"
#include "util_i.h"
#include "util.h"
#include "util_genes.h"

#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[1000];
	int i, j;
	int num_genes = 0;
	int num_exons = 0;
	struct g_list *genes;
	struct exons_list *exons;
	int b, e;
	char name[100], scf_name[100];
	int cur_exons_count;

	if( argc != 2 ) {
		printf("sort_exons exons_file\n");
		return EXIT_FAILURE;
	}

	strcpy(name, "");
	strcpy(scf_name, "");
	f = fopen(argv[1], "r");

	while(fgets(buf, 1000, f))
	{
		if( buf[0] == '#' ) {}	
		else if((buf[0] == '>') || (buf[0] == '<')) num_genes++;
		else num_exons++;

	}

	if( num_genes > 0 ) {
		genes = (struct g_list *) ckalloc(sizeof(struct g_list) * num_genes);
	}
	else {
		genes = (struct g_list *) ckalloc(sizeof(struct g_list));
	}

	if( num_exons > 0 ) {
		exons = (struct exons_list *) ckalloc(sizeof(struct exons_list) * num_exons);
	}
	else {
		exons = (struct exons_list *) ckalloc(sizeof(struct exons_list));
	}
	
	fseek(f, 0, SEEK_SET);

	i = -1;
	j = 0;
	while(fgets(buf, 1000, f))
	{
		if( buf[0] == '#' ) {}
		else if( (buf[0] == '>') || (buf[0] == '<') )
		{
			if( i >= 0 ) {
				genes[i].exonCount = cur_exons_count;
				genes[i].exEnd = j-1;
			}
			i++;
			cur_exons_count = 0;
			if(buf[0] == '>') genes[i].strand = '+';
			else if(buf[0] == '<' ) genes[i].strand = '-';
			else fatalf("unexpected strand %c\n", buf[0]);

			if( sscanf(buf, "%*s %d %d %s %s %*s", &b, &e, name, scf_name) == 4 ) {
				strcpy(genes[i].sname, scf_name);
			}
			else if( sscanf(buf,  "%*s %d %d %s %*s", &b, &e, name) != 3 ) {
				printf("wrong format in %s\n", buf);	
			}
			else strcpy(genes[i].sname, "");

			genes[i].gid = i;
			genes[i].txStart = b;
			genes[i].txEnd = e;
			strcpy(genes[i].gname, name);
		}
		else {
			sscanf(buf, "%d %d", &b, &e);
			if( cur_exons_count == 0 ) genes[i].exStart = j;
			exons[j].fid = i;
			exons[j].reg = assign_I(b, e);
			cur_exons_count++;
			j++;
		}
	}
	genes[i].exonCount = cur_exons_count;
	genes[i].exEnd = j-1;

	quick_sort_inc_genes(genes, 0, num_genes-1, POS_BASE);
	i = 0;
	while( i < num_genes ) {
		j = 0;
		while( ((i+j) < num_genes) && (genes[i].txStart == genes[i+j].txStart )) j++;
		quick_sort_dec_genes(genes, i, i+j-1, LEN_BASE);
		i = i+j;
	}

	for( i = 0; i < num_genes; i++ ) {
		if( genes[i].txStart < 0 ) {} 
		else {
			if( genes[i].strand == '+' ) {
				if( strcmp(genes[i].sname, "") == 0 ) {
					printf("> %d %d %s\n", genes[i].txStart, genes[i].txEnd, genes[i].gname);
				}
				else {
					printf("> %d %d %s %s\n", genes[i].txStart, genes[i].txEnd, genes[i].gname, genes[i].sname);

				}
			}
			else if( genes[i].strand == '-' ) {
				if( strcmp(genes[i].sname, "") == 0 ) {
					printf("< %d %d %s (complement)\n", genes[i].txStart, genes[i].txEnd, genes[i].gname);
				}
				else {
					printf("< %d %d %s %s (complement)\n", genes[i].txStart, genes[i].txEnd, genes[i].gname, genes[i].sname);
				}
			}
			else fatalf("unexpected strand %c\n", genes[i].strand);

			for( j = genes[i].exStart; j <= genes[i].exEnd; j++ ) {
				printf("%d %d\n", exons[j].reg.lower, exons[j].reg.upper);
			}
		}
	}
	fclose(f);

	free(genes);
	free(exons);

	return EXIT_SUCCESS;
}
