#include "sort_genes.h"
#include "util_genes.h"
#include "util.h"
#include "util_i.h"

extern int debug_mode;

void quick_sort_dec_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].txStart));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].txEnd - a[(lo+hi)/2].txStart);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].txStart))>x) i++; 
			while (((float)(a[j].txStart))<x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].txEnd-a[i].txStart)>x) i++; 
			while (abs(a[j].txEnd-a[j].txStart)<x) j--;
		}

		if (i<=j)
		{
			h = assign_genes(a[i]);
			strcpy(h.gname, a[i].gname);
			strcpy(h.sname, a[i].sname);
			a[i] = assign_genes(a[j]);
			strcpy(a[i].gname, a[j].gname);
			strcpy(a[i].sname, a[j].sname);
			a[j] = assign_genes(h);
			strcpy(a[j].gname, h.gname);
			strcpy(a[j].sname, h.sname);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_genes(a, lo, j, mode);
	if (i < hi) quick_sort_dec_genes(a, i, hi, mode);
}

void quick_sort_inc_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].txStart));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].txEnd - a[(lo+hi)/2].txStart);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].txStart))<x) i++; 
			while (((float)(a[j].txStart))>x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].txEnd-a[i].txStart)<x) i++; 
			while (abs(a[j].txEnd-a[j].txStart)>x) j--;
		}

		if (i<=j)
		{
			h = assign_genes(a[i]);
			strcpy(h.gname, a[i].gname);
			strcpy(h.sname, a[i].sname);
			a[i] = assign_genes(a[j]);
			strcpy(a[i].gname, a[j].gname);
			strcpy(a[i].sname, a[j].sname);
			a[j] = assign_genes(h);
			strcpy(a[j].gname, h.gname);
			strcpy(a[j].sname, h.sname);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_genes(a, lo, j, mode);
	if (i < hi) quick_sort_inc_genes(a, i, hi, mode);
}

int quick_search_close_genes(struct g_list *sorted, int i, int j, int query)
{
	int mid;
	int val;
	int res;

	mid = (i+j)/2;
	val = sorted[mid].txStart;

	if(val > query) {
		if( j <= (mid+1) ) return (mid+1);	
		else res = quick_search_close_genes(sorted, mid+1, j, query);
	} 
	else if(val < query) {
		if( i >= (mid-1) ) return i;
		else res = quick_search_close_genes(sorted, i, mid-1, query);
	}
	else return mid;

	return(res);
}

struct g_list assign_genes(struct g_list a)
{
  struct g_list res;

  res.gid = a.gid;
  res.strand = a.strand;
  res.txStart = a.txStart;
  res.txEnd = a.txEnd;
  res.exonCount = a.exonCount;
  res.exStart = a.exStart;
  res.exEnd = a.exEnd;
	strcpy(res.gname, a.gname);

  return(res);
}
