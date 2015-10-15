#include "bed.h"
#include "util_i.h"

void init_bed(struct bed *list, int num_rows)
{
	int i = 0;

	for( i = 0; i < num_rows; i++ ) {
		strcpy(list[i].chr, ".");	
		list[i].reg = assign_I(0, 1);
		list[i].qual = (float) 0;
		strcpy(list[i].feature, ".");	
		strcpy(list[i].ref, ".");	
		strcpy(list[i].alt, ".");	
		strcpy(list[i].pass, ".");	
		strcpy(list[i].info, ".");	
	}
}
