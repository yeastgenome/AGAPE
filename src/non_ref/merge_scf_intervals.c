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
	int b = 0, e = 1;
	int count = 0;
	int num_match_regions = 0, num_merged = 0;
	struct I *match_regions, *merged_regions;
	char scaf_name[MAX_NAME], cur_name[MAX_NAME];
	char buf[MAX_NAME];
	char strain_name[MAX_NAME];
	bool is_pid = false;
	int *pid;
	int cur_pid = 0;
	bool size_print = false;
	bool is_strain_name_given = false;
	int size = 0, old_size = 0;

	debug_mode = FALSE;
	if( argc == 3 ) {
		if( strcmp(argv[2], "Size_Print") == 0 ) {
			size_print = true;
		}
		else {
			strcpy(strain_name, argv[2]);
			is_strain_name_given = true;
		}
	}
	else if( argc != 2 ) {
		fatal("args: maf scaf_name\n");
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
		match_regions = (struct I *) ckalloc(count * (sizeof(struct I)) );
		merged_regions = (struct I *) ckalloc(count * (sizeof(struct I)) );
		pid = (int *) ckalloc(count * (sizeof(int)) );
		initialize_I_list(match_regions, count);
		initialize_I_list(merged_regions, count);

		fseek(f, 0, SEEK_SET);
		while( fgets(buf, MAX_NAME, f) ) {
			old_size = size;
			if( sscanf(buf, "%s %d %d %d", cur_name, &b, &e, &cur_pid) == 4 ) {
				if( size_print == true ) {
					size = cur_pid;
				}
				else {
					is_pid = true;
				}
			}
			else if( sscanf(buf, "%s %d %d", cur_name, &b, &e) != 3 ) {
				fatalf("wrong line : %s", buf);
			}

			if( (k == 0) || (strcmp( cur_name, scaf_name ) == 0) ) {
				if( is_pid == true ) pid[i] = cur_pid;
				match_regions[i] = assign_I(b, e);
				i++;
			}
			else {
				if( i > 1 ) {		
					num_match_regions = i;
					if( is_pid == true ) {
						num_merged = sort_merge_intervals_and_pid(match_regions, num_match_regions, merged_regions, pid);
					}
					else {
						num_merged = sort_merge_intervals(match_regions, num_match_regions, merged_regions);
					}
					for( j = 0; j < num_merged; j++ ) {
						if( is_strain_name_given == true ) {
							if( is_pid == true ) {
								printf("%s %d %d %s %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, strain_name, pid[j]);
							}
							else {
								printf("%s %d %d %s\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, strain_name);
							}
						}
						else {
							if( size_print == true ) {
								printf("%s %d %d %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, old_size);
							}
							else {
								printf("%s %d %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper);
							}
						}
					}
				}
				else {
					if( is_strain_name_given == 3 ) {
						if( is_pid == true ) {
							printf("%s %d %d %s %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper, strain_name, pid[0]);
						}
						else {
							printf("%s %d %d %s\n", scaf_name, match_regions[0].lower, match_regions[0].upper, strain_name);
						}
					}
					else {
						if( size_print == true ) {
							printf("%s %d %d %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper, old_size);
						}
						else {
							printf("%s %d %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper);
						}
					}
				}

				i = 0;
				match_regions[i] = assign_I(b, e);
				if( is_pid == true ) {
					pid[i] = cur_pid;
				}
				i++;
			}	
			strcpy(scaf_name, cur_name);
			k++;
		}
	}

	if( i > 1 ) {		
		num_match_regions = i;
		if( is_pid == true ) {
			num_merged = sort_merge_intervals_and_pid(match_regions, num_match_regions, merged_regions, pid);
		}
		else {
			num_merged = sort_merge_intervals(match_regions, num_match_regions, merged_regions);
		}

		for( j = 0; j < num_merged; j++ ) {
			if( is_strain_name_given == true ) {
				if( is_pid == true ) {
					printf("%s %d %d %s %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, strain_name, pid[j]);
				}
				else {
					printf("%s %d %d %s\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, strain_name);
				}
			}
			else {
				if( size_print == true ) {
					printf("%s %d %d %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper, size);
				}
				else {
					printf("%s %d %d\n", scaf_name, merged_regions[j].lower, merged_regions[j].upper);
				}
				
			}
		}
	}
	else {
		if( is_strain_name_given == true ) {
			if( is_pid == true ) {
				printf("%s %d %d %s %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper, strain_name, pid[0]);
			}
			else {
				printf("%s %d %d %s\n", scaf_name, match_regions[0].lower, match_regions[0].upper, strain_name);
			}
		}
		else {
			if( size_print == true ) {
				printf("%s %d %d %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper, size);
			}
			else {
				printf("%s %d %d\n", scaf_name, match_regions[0].lower, match_regions[0].upper);
			}
		}
	}	

	if( count > 0 ) {
		free(match_regions);
		free(merged_regions);
	}
	return EXIT_SUCCESS;
}
