// interval: [a,b]
#include "main.h"
#include "util.h"
#include "input_blastn.h"

#define INTERVAL_LEN 100
#define LEN_RATIO 0.2
#define LEN_FRACTION 0.8
#define HOMOLOGOUS_MATCH_PID 80

int debug_mode;

int main(int argc, char **argv)
{
	struct sp_list **homologs;
	int num_homologs = 0;
	int *max_hits;
	FILE *f;
	int i = 0, j = 0;
	int count = 0;
	bool is_print = true;

	if( argc != 2 ) {
		fatal("args: blastn_out\n");
	}

	max_hits = (int *) ckalloc(sizeof(int));
	*max_hits = 1;
	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[1]);    
	}
	else {
		count = count_genes_blastn(f, max_hits);
	}

	j = 0;
	if( count > 0 ) {
		homologs = (struct sp_list **) ckalloc(count * (sizeof(struct sp_list *)) );
		for( i = 0; i < count; i++ ) {
			homologs[i] = (struct sp_list *) ckalloc((*max_hits) * (sizeof(struct sp_list)));
			for( j = 0; j < (*max_hits); j++ ) {
				homologs[i][j].id = 0;
				strcpy(homologs[i][j].name, "");
			}
		}
 		num_homologs = input_genes_blastn(f, homologs, HOMOLOGOUS_MATCH_PID, is_print);
	}	

	free(max_hits);
	if( count > 0 ) {
		for( i = 0; i < count; i++ ) free(homologs[i]);
		free(homologs);
	}
	fclose(f);
	return EXIT_SUCCESS;
}
