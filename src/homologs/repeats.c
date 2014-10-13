#include "main.h"
#include "util_exons.h"
#include "util_i.h"
#include "regions.h"
#include "repeats.h"

bool is_repeats(struct exons_list *exons, int num_exons, char *name, int from, int to) // exons assume to be already sorted by genomic positions
{
	int i = 0;
	int mid = 0;
	bool res = false;
	struct I reg;

	if( to > from ) {
		reg = assign_I(from, to);
	}
	else {
		fatalf("unexpected interval: %d-%d\n", from, to);
	}

	mid = quick_search_close_exons(exons, 0, num_exons-1, from);
	i = mid;
	while( (res == false) && (i < num_exons) && (exons[i].reg.lower <= to)) {
		if( width(reg) <= SHORT_LEN_TH ) {
// (strcmp(reg, exons[i].reg) == 0) && (almost_subset(reg, exons[i].reg) == true) ) {
//			printf("%d-%d too short\n", from, to);
			res = true;	
		}
		else if( (strcmp(name, exons[i].chr) == 0) && subset(reg, exons[i].reg) ) {
//			printf("%d-%d belongs to %s %d-%d\n", from, to, exons[i].chr, exons[i].reg.lower, exons[i].reg.upper);
			res = true;
		}
		i++;
	}

	return(res);
}
