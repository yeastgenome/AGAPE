#include "main.h"
#include "util_i.h"
#include "util.h"

bool proper_overlap(struct I reg1, struct I reg2)
{
  if(overlap(reg1, reg2))
  {
    if(width(intersect(reg1, reg2)) == 0) return false;
    else return true;
  }
  else return false;
}

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

	if( lo > hi ) fatalf("empty interval in (%d,%d)", lo, hi);
	else {
		res.lower = lo;
		res.upper = hi;
		return(res);
	}
}

bool near_equal(struct I a, struct I b, int th)
{
	bool res = false;

	if((abs(a.lower-b.lower) <= th) && (abs(a.upper-b.upper) <= th) )
	{
		res = true;
	}
	else res = false;

	return(res);
}

bool mostly_overlap(struct I a, struct I b, float cut_ratio)
{
	bool res = false;
	int min_width = 0;
	struct I reg, min_reg;

	reg = assign_I(0, 1);
	min_reg = assign_I(0, 1);

	if( width(a) > width(b) ) {
		min_width = width(b);
		min_reg = assign_I(b.lower, b.upper);
	}
	else {
		min_width = width(a);
		min_reg = assign_I(a.lower, a.upper);
	}

	reg = intersect(a, b);

	if( (float)width(reg)/(float)width(min_reg) >= cut_ratio )
	{
		res = true;
	}
	else res = false;

	return(res);
}
