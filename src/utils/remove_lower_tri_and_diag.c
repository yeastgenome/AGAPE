/* remove_lower_tri_and_diagonal - remove redundant and trivial alignments */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maf.h"
#include "util.h"
#include "mz_scores.h"

#define MAX_LEN 2000

int main(int argc, char** argv) {
	struct mafFile* mf;
	struct mafAli* ali;
	struct mafComp* mc;
	FILE *f;	
	char buf[MAX_LEN];

	if ( argc != 2) {
		printf("remove_lower_tri_and_diag maf-file\n");
		return EXIT_FAILURE;
	}
	
//	init_scores70();
//	mafWriteStart(stdout, 0);
	
	f = ckopen(argv[1], "r");
	while(fgets(buf, MAX_LEN, f)) {
		if( buf[0] == '#' ) {
			printf("%s", buf);
		}
		else break;
	}	

	fseek(f, 0, SEEK_SET);
	fclose(f);
	
	mf = mafOpen(argv[1], 0);
	while((ali = mafNext(mf)) != NULL) {
		mc = ali->components;

		if(mc->next->strand == '+' && mc->start > mc->next->start)
			continue;
		else if(mc->next->strand == '-' && mc->start > (mc->next->srcSize - mc->next->start - mc->next->size))
			continue;
		else if((strcmp(mc->name, mc->next->name) == 0) && (mc->next->strand == '+') && (mc->start == mc->next->start) && (mc->size == mc->next->size))
			continue;
		else
			mafWrite(stdout, ali);
	}

	mafFileFree(&mf);
	
	mafWriteEnd(stdout);

	return EXIT_SUCCESS;
}

