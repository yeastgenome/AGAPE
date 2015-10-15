#include "main.h"
#include "util_i.h"
#include "util.h"

bool overlap(struct I a, struct I b)
{
	if( (a.lower >= b.lower) && (a.lower <= b.upper) )
	{
		return(true);
	}
	else if( (a.upper >= b.lower) && (a.upper <= b.upper) )
	{
		return(true);
	}
	else if( (b.lower >= a.lower) && (b.lower <= a.upper))
	{
		return(true);
	}
	else if( (b.upper >= a.lower) && (b.upper <= a.upper))
	{
		return(true);
	}
	else return(false);
}

int width(struct I temp)
{
	if( temp.lower > temp.upper ) {
		fatalf("empty interval in (%d,%d)", temp.lower, temp.upper);
	}
	else return(temp.upper - temp.lower);
}

struct I assign_I(int lower, int upper)
{
	struct I res;

	res.lower = lower;
	res.upper = upper;

	return(res);
}

bool subset(struct I a, struct I b)
{
	if( (a.lower >= b.lower) && (a.upper <= b.upper) ) return(true);
	else return(false);
}

bool proper_subset(struct I a, struct I b)
{
	if( (a.lower > b.lower) && (a.upper < b.upper) ) return(true);
	else return(false);
}

bool equal(struct I a, struct I b)
{
	if( (a.lower == b.lower) && (a.upper == b.upper) ) return(true);
	else return(false);
}

bool in(int r, struct I reg)
{
	if( (r <= reg.upper) && (r >= reg.lower) ) return(true);
	else return(false);
}

struct I intersect(struct I a, struct I b)
{
	struct I res;
	int lo, hi;

	if( a.lower >= b.lower ) lo = a.lower;
	else lo = b.lower;

	if( a.upper <= b.upper ) hi = a.upper;
	else hi = b.upper;

	if( lo > hi ) {
		fatalf("empty interval in (%d,%d)", lo, hi);
	}
	else {
		res.lower = lo;
		res.upper = hi;
		return(res);
	}
}

