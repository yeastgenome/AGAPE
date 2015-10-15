// interval: [a,b]
#include "main.h"
#include "read_maf.h"
#include "util_I_gen.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"
#include "regions.h"

#define INTERVAL_LEN 300

#define REF_MATCH_PID 70

bool is_in_name_list(char *name, char **list, int num_items);

bool debug_mode;
char S[BIG], T[BIG];

int main(int argc, char **argv)
{
	struct DotList *algns;
	int *num_algns;
	int *size1, *size2;
	FILE *f;
	int i = 0, j = 0;
	int b = 0, e = 0;
	int count = 0;
	char species[MAX_NAME], species2[MAX_NAME];
	bool is_end = false;
	int len_cutoff = 0;
	int pid_cutoff = 0;
	int *match_pid;
	struct I *temp_intervals;
	struct I *new_intervals;
	int num_intervals = 0;
	int contig_len = 0;
	char contig_name[MAX_NAME];

	debug_mode = false;
	pid_cutoff = REF_MATCH_PID;
	len_cutoff = INTERVAL_LEN;
	if( argc != 4 ) {
		fatal("args: maf contig_name contig_length\n");
	}
	else {
		strcpy(contig_name, argv[2]);
		contig_len = atoi(argv[3]);
	}

	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	strcpy(S, "");
	strcpy(T, "");

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
		match_pid = (int *) ckalloc(count * sizeof(int));
		for( i = 0; i < count; i++ ) {
			match_pid[i] = 0;
		}
		temp_intervals = (struct I *) ckalloc(count * sizeof(struct I));
		new_intervals = (struct I *) ckalloc(count * sizeof(struct I));
		initialize_I_list(temp_intervals, count);
		initialize_I_list(new_intervals, count);

		for( i = 0; i < (*num_algns); i++ ) {
			if( (algns[i].identity >= pid_cutoff) && (strcmp(algns[i].name2, "chrM") != 0) ) {
				temp_intervals[j] = assign_I(algns[i].x.lower, algns[i].x.upper);
				match_pid[j] = algns[i].identity;
				j++;
			}
		}

		num_intervals = sort_merge_intervals_and_pid(temp_intervals, j, new_intervals, match_pid);
		for( j = 0; j < num_intervals; j++ ) {
			if( j == 0 ) b = 1;
			else {
				b = new_intervals[j-1].upper + 1;
				e = new_intervals[j].lower - 1;
			}
			if( (e-b) > len_cutoff ) {
				printf("%s %d %d\n", contig_name, b, e);
			}
		}
		b = new_intervals[num_intervals-1].upper + 1;
		e = contig_len;
		if( (e-b) > len_cutoff ) {
			printf("%s %d %d\n", contig_name, b, e);
		}
	}	

	free(num_algns);
	if( count > 0 ) {
		free(algns);
		free(match_pid);
		free(temp_intervals);
		free(new_intervals);
	}
	free(size1);
	free(size2);
	return EXIT_SUCCESS;
}

bool is_in_name_list(char *name, char **list, int num_items)
{
	int i = 0;
	bool res = false;

	while( ( i < num_items) && (res == false) ) {
		if( strcmp(name, list[i]) == 0 ) {
			res = true;
		}
		i++;
	}
	return(res);
}
