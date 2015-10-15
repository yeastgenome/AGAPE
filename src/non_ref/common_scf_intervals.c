// interval: [a,b]
#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_I_gen.h"
#include "util_gen.h"
#include "regions.h"

#define SHORT_LEN 50
#define LONG_LEN 100
#define INTERVAL_LEN 300

#define HIGH_PID_IN_LONG 75
#define HIGH_PID_IN_SHORT 90

int debug_mode;

int main(int argc, char **argv)
{
	FILE *f;
	int i = 0, j = 0, k = 0;
	int count = 0;
	int num_match_regions1 = 0, num_merged = 0;
	int num_match_regions2 = 0;
	struct scf_I *match_regions1;
	struct scf_I *match_regions2;
	struct I *merged_regions;
	char scaf_name[MAX_NAME], cur_name[MAX_NAME];
	char buf[MAX_NAME];
	struct I *cur_regions;

	debug_mode = FALSE;
	if( argc == 4 ) {
		debug_mode = TRUE;
	}
	else if( argc != 3 ) {
		fatal("args: gff1 gff2\n");
	}

	strcpy(buf, "");
	strcpy(scaf_name, "");
	strcpy(cur_name, "");

	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[1]);    
	}
	else {
		while(fgets(buf, MAX_NAME, f)) count++;
	}

	if( count > 0 ) {
		match_regions1 = (struct scf_I *) ckalloc(count * (sizeof(struct scf_I)) );
		initialize_scf_I_list(match_regions1, count);
		num_match_regions1 = input_scf_I_list(f, match_regions1, count);
	}
	fclose(f);

	count = 0;
	if( (f = fopen(argv[2], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[2]);    
	}
	else {
		while(fgets(buf, MAX_NAME, f)) count++;
	}

	if( count > 0 ) {
		match_regions2 = (struct scf_I *) ckalloc(count * (sizeof(struct scf_I)) );
		initialize_scf_I_list(match_regions2, count);
		num_match_regions2 = input_scf_I_list(f, match_regions2, count);
	}
	fclose(f);

	count = num_match_regions1 + num_match_regions2;
	if( count > 0 ) {
		cur_regions = (struct I * ) ckalloc( count * sizeof(struct I));
		merged_regions = (struct I * ) ckalloc( count * sizeof(struct I));
		initialize_I_list(cur_regions, count);
		initialize_I_list(merged_regions, count);
	}
	i = 0;
	strcpy(scaf_name, "");
	while( i < num_match_regions1 ) {	
		strcpy(cur_name, match_regions1[i].name);
		if( i == 0 ) {
			strcpy(scaf_name, cur_name);
			cur_regions[j] = assign_I(match_regions1[i].region.lower, match_regions1[i].region.upper);
			j++;
		}
		else if( strcmp(scaf_name, cur_name) == 0 ) {
			cur_regions[j] = assign_I(match_regions1[i].region.lower, match_regions1[i].region.upper);
			j++;
		}
		else { // scaf_name != cur_name
			count = find_match_regions(cur_regions, j, match_regions2, num_match_regions2, scaf_name);
			if( count > 0 ) {
				count = count + j;
				num_merged = sort_merge_intervals(cur_regions, count, merged_regions);	
				for( k = 0; k < num_merged; k++ ) {
					printf("%s %d %d\n", scaf_name, merged_regions[k].lower, merged_regions[k].upper);
				}
			}
			j = 0;
			strcpy(scaf_name, cur_name);
			cur_regions[j] = assign_I(match_regions1[i].region.lower, match_regions1[i].region.upper);
			j++;
		}
		i++;	
	}

	if( count > 0 ) {
		free(match_regions1);
		free(match_regions2);
		free(cur_regions);
		free(merged_regions);
	}
	return EXIT_SUCCESS;
}
