// interval: [a,b]
#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "util_I_gen.h"
#include "util_gen.h"
#include "regions.h"

int debug_mode;
char S[BIG], T[BIG];

int main(int argc, char **argv)
{
	struct DotList *algns;
	int *num_algns;
	int *size1, *size2;
	FILE *f;
	int i = 0, j = 0;
	int count = 0;
	int num_match_regions = 0, num_merged = 0;
	struct I *match_regions, *merged_regions;
	char scf_name[MAX_NAME];
	char chr_name[MAX_NAME];
	char species[MAX_NAME], species2[MAX_NAME];
	bool is_end = false;
	int total_len = 0;
	float ratio = (float) 0;
	int b = 0, e = 1;

	debug_mode = FALSE;
	if( argc == 7 ) {
		debug_mode = TRUE;
	}
	else if( argc != 6 ) {
		fatal("args: maf chr b e fraction\n");
	}

	strcpy(chr_name, "");
	strcpy(chr_name, argv[2]);
	b = atoi(argv[3]);
	e = atoi(argv[4]);
	ratio = atof(argv[5]);
	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	strcpy(S, "");
	strcpy(T, "");
	strcpy(scf_name, "");

	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[1]);    
	}
	else {
		while((is_end == false) && fgets(S, BIG, f)) {
			if( S[0] == '#' ) {
				while( (is_end == false) && (S[0] == '#') ) {
					if( fgets(S, BIG, f) == NULL ) is_end = true;
					if( strncmp(S, "##maf", 5) == 0 ) {
						count = 0;
					}
				}
				count = 0;	
			}

 		 	if( S[0] == 'a' ) {
				count++;
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", species) != 1) || (sscanf(T, "%*s %s %*s", species2) != 1)) {}
			}
		}
	}
	fclose(f);

	j = 0;
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * (sizeof(struct DotList)) );
 		read_maf(argv[1], G_MODE, algns, num_algns, size1, size2);
		strcpy(scf_name, species);
		match_regions = (struct I *) ckalloc(count * sizeof(struct I));
		initialize_I_list(match_regions, count);

		for( i = 0; i < (*num_algns); i++ ) {
			if( strcmp(algns[i].name2, chr_name) == 0 ) {
				match_regions[j] = assign_I(algns[i].x.lower, algns[i].x.upper);	
				j++;
			}
		}
		
		num_match_regions = j;
		if( num_match_regions > 0 ) {
			merged_regions = (struct I *) ckalloc(num_match_regions * sizeof(struct I));
			initialize_I_list(merged_regions, num_match_regions);
			num_merged = sort_merge_intervals(match_regions, num_match_regions, merged_regions);
			free(match_regions);
		}
	}	

	total_len = 0;
	for( i = 0; i < num_merged; i++ ) {
		total_len = total_len + width(merged_regions[i]);
	}

	if( ((float)total_len/(float)(*size1)) < ratio ) 
	{
		printf("%s %d %d\n", scf_name, b, e);
	} 

	free(num_algns);
	if( count > 0 ) {
		free(algns);
	}
	free(size1);
	free(size2);
	if( num_match_regions > 0 ) {
		free(merged_regions);
	}
	return EXIT_SUCCESS;
}
