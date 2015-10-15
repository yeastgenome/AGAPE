// interval: [a,b]
#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"
#include "regions.h"

#define INTERVAL_LEN 300

#define IDENTICAL_MATCH_PID 99

bool debug_mode;
char S[BIG], T[BIG];

int main(int argc, char **argv)
{
	struct DotList *algns;
	int *num_algns;
	int *size1, *size2;
	FILE *f;
	int i = 0, j = 0;
	int count = 0;
	char species[MAX_NAME], species2[MAX_NAME];
	bool is_end = false;
	int len_cutoff = 0;
	int pid_cutoff = 0;
	int Size_Print;

	debug_mode = false;
	pid_cutoff = IDENTICAL_MATCH_PID;
	Size_Print = FALSE;
	len_cutoff = INTERVAL_LEN;
	if( argc == 5 ) {
		Size_Print = TRUE;
		len_cutoff = (int)atof(argv[2]);
		pid_cutoff = atoi(argv[3]);
	}
	else if( argc == 4 ) {
		len_cutoff = atoi(argv[2]);
		pid_cutoff = atoi(argv[3]);
	}
	else if( argc == 3 ) {
		pid_cutoff = atoi(argv[2]);
	}
	else if( argc != 2 ) {
		fatal("args: maf\n");
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

		for( i = 0; i < (*num_algns); i++ ) {
			if( width(algns[i].x) > len_cutoff ) {
				if( algns[i].identity >= pid_cutoff ) {
					if( Size_Print == TRUE ) {
						printf("%s %d %d %d\n", algns[i].name1, algns[i].x.lower, algns[i].x.upper, algns[i].len1);
					}
					else {
						printf("%s %d %d\n", algns[i].name1, algns[i].x.lower, algns[i].x.upper);
					}
				}
			}
		}
	}	

	free(num_algns);
	if( count > 0 ) {
		free(algns);
	}
	free(size1);
	free(size2);
	return EXIT_SUCCESS;
}

