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
	int i = 0;
	int count = 0;
	int num_match_regions = 0;
	struct orf_I * match_regions;
	char scaf_name[MAX_NAME], cur_name[MAX_NAME];
	char buf[MAX_NAME];
	struct I reg;
	int b = 0, e = 0;

	reg = assign_I(0, 1);

	debug_mode = FALSE;
	if( argc == 4 ) {
		debug_mode = TRUE;
	}
	else if( argc != 3 ) {
		fatal("args: intervals1 intervals2\n");
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
		match_regions = (struct orf_I *) ckalloc(count * (sizeof(struct orf_I)) );
		initialize_orf_I_list(match_regions, count);
		num_match_regions = input_orf_I_list(f, match_regions, count);
	}
	fclose(f);

	count = 0;
	if( (f = fopen(argv[2], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[2]);    
	}
	else {
		while(fgets(buf, MAX_NAME, f)) {
			if( buf[0] == '>' ) {
				printf("%s", buf);
			}
			else {
				if( sscanf(buf, "%s %d %d %s %*s", scaf_name, &b, &e, cur_name) != 4 ) {
					fatalf("wrong interval line: %s", buf);    
				}
				else {
					i = 0;
					reg = assign_I(b, e);
					while( i < num_match_regions ) {
						if( strcmp(cur_name, match_regions[i].strain_name) == 0 ) {
							if( strcmp(scaf_name, match_regions[i].name) == 0 ) {
								if( proper_overlap(reg, match_regions[i].region) == true ) {
									printf("%s %d %d %s\n", match_regions[i].name, match_regions[i].region.lower, match_regions[i].region.upper, match_regions[i].strain_name);
								}							
							}
						}	
						i++;
					}
				}
			}
		}
	}

	if( count > 0 ) {
		free(match_regions);
	}
	return EXIT_SUCCESS;
}
