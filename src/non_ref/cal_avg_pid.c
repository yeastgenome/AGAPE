// interval: [a,b]
#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "util_I_gen.h"
#include "util_gen.h"
#include "regions.h"

#define SHORT_LEN 50
#define LONG_LEN 100

#define HIGH_PID_IN_LONG 75
#define HIGH_PID_IN_SHORT 90

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
	char species[MAX_NAME], species2[MAX_NAME];

	debug_mode = FALSE;
	if( argc == 3 ) {
		debug_mode = TRUE;
	}
	else if( argc != 2 ) {
		fatal("args: maf scaf_name\n");
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
		while(fgets(S, BIG, f)) {
			if( S[0] == '#' ) {
				while( S[0] == '#' ) {
					fgets(S, BIG, f);
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
			if( width(algns[i].x) > LONG_LEN ) {
				printf("%d %d\n", width(algns[i].x), algns[i].identity);
			}
		}
	}
		
	free(num_algns);
	free(size1);
	free(size2);
	if( count > 0 ) {
		free(algns);
	}
	return EXIT_SUCCESS;
}
