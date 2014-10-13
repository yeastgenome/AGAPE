#include "util_sort.h"
#include "util.h"

void quick_sort_dec(struct slist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct slist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while ((i <= hi) && (a[i].val>x)) i++; 
		while ((j >= lo) && (a[j].val<x)) j--;
		if (i<=j)
		{
			h = assign_slist(a[i]);
			a[i] = assign_slist(a[j]);
			a[j] = assign_slist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec(a, lo, j);
	if (i < hi) quick_sort_dec(a, i, hi);
}

void quick_sort_inc(struct slist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct slist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while ((i <= hi) && (a[i].val<x)) i++; 
		while ((j >= lo) && (a[j].val>x)) j--;
		if (i<=j)
		{
			h = assign_slist(a[i]);
			a[i] = assign_slist(a[j]);
			a[j] = assign_slist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc(a, lo, j);
	if (i < hi) quick_sort_inc(a, i, hi);
}

struct slist assign_slist(struct slist a)
{
	struct slist res;

	res.id = a.id;
	res.val = a.val;

	return(res);
}
