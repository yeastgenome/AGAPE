#include "util_exons.h"
#include "util_i.h"
#include "util.h"
#include "regions.h"

void init_exons(struct exons_list *exons, int from, int to)
{
	int i = 0;

	for( i = from; i <= to; i++ ) {
		exons[i].reg = assign_I(0, 1);
		exons[i].dir = '+'; // '+' or '-'
		strcpy(exons[i].chr, ""); // for repeats, a chromosome name is saved
		strcpy(exons[i].name, ""); 
	}
}

void assign_gff_exons(FILE *f, struct exons_list *exons, int num_exons)
{
	assign_gff_exons_mode(f, exons, num_exons, "", ALL_MODE);
}

void assign_gff_exons_chr(FILE *f, struct exons_list *exons, int num_exons, char *chr)
{
	assign_gff_exons_mode(f, exons, num_exons, chr, CHR_MODE);
}

void assign_gff_exons_mode(FILE *f, struct exons_list *exons, int num_exons, char *chr, int mode)
{
	char *buf;
	int i = 0;
	int b = 0, e = 0;
	char sign = '+';
	char name[LEN_NAME];

	buf = (char *) ckalloc(NUM_CHARS * sizeof(char));
	while(fgets(buf, NUM_CHARS, f)){
		if( sscanf(buf, "%s %*s %*s %d %d %*s %c %*s", name, &b, &e, &sign) != 4) {
			fatalf("format error in %s\n", buf);
		}
		else {
			if( (mode == ALL_MODE) || ( (mode == CHR_MODE) && (strcmp(name, chr) == 0) ) ) 
			{
				strcpy(exons[i].chr, name);
				exons[i].dir = sign;
				exons[i].reg = assign_I(b, e);
			}
		}
		i++;
		if( i > num_exons) {
			fatalf("exon counting error: exceed %d\n", num_exons);
		}
	}
}

void quick_sort_dec_exons(struct exons_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct exons_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].reg.lower));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].reg.upper - a[(lo+hi)/2].reg.lower);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].reg.lower))>x) i++; 
			while (((float)(a[j].reg.lower))<x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].reg.upper-a[i].reg.lower)>x) i++; 
			while (abs(a[j].reg.upper-a[j].reg.lower)<x) j--;
		}

		if (i<=j)
		{
			h = assign_exons(a[i]);
			a[i] = assign_exons(a[j]);
			a[j] = assign_exons(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_exons(a, lo, j, mode);
	if (i < hi) quick_sort_dec_exons(a, i, hi, mode);
}

void quick_sort_inc_exons(struct exons_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct exons_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].reg.lower));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].reg.upper - a[(lo+hi)/2].reg.lower);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].reg.lower))<x) i++; 
			while (((float)(a[j].reg.lower))>x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].reg.upper-a[i].reg.lower)<x) i++; 
			while (abs(a[j].reg.upper-a[j].reg.lower)>x) j--;
		}

		if (i<=j)
		{
			h = assign_exons(a[i]);
			a[i] = assign_exons(a[j]);
			a[j] = assign_exons(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_exons(a, lo, j, mode);
	if (i < hi) quick_sort_inc_exons(a, i, hi, mode);
}

void selection_sort_exons(struct exons_list *a, int num)
{
	int min = 0, i = 0, j = 0;
	struct exons_list temp;

	for( i = 0; i < num; i++ ) {
		min = i;
		for( j = i+1; j < num; j++ ) {
			if( strcmp(a[j].chr, a[min].chr) < 0 ) {
				min = j;	
			}
			else if( strcmp(a[j].chr, a[min].chr) == 0 ) {
				if( a[j].reg.lower < a[min].reg.lower ) {
					min = j;
				}
				else if( a[j].reg.lower == a[min].reg.lower ) {
					if( width(a[j].reg) > width(a[min].reg) ) {
						min = j;
					}
				}
			}
		}	
		temp = assign_exons(a[min]);
		a[min] = assign_exons(a[i]);
		a[i] = assign_exons(temp);
	}
}

int quick_search_close_exons(struct exons_list *sorted, int i, int j, int query)
{
	int mid;
	int val;
	int res;

	mid = (i+j)/2;
	val = sorted[mid].reg.lower;

	if(val > query) {
		if( j <= (mid+1) ) return (mid+1);	
		else res = quick_search_close_exons(sorted, mid+1, j, query);
	} 
	else if(val < query) {
		if( i >= (mid-1) ) return i;
		else res = quick_search_close_exons(sorted, i, mid-1, query);
	}
	else return mid;

	return(res);
}

struct exons_list assign_exons(struct exons_list a)
{
  struct exons_list res;

  res.dir = a.dir;
  res.reg.lower = a.reg.lower;
  res.reg.upper = a.reg.upper;
	strcpy(res.chr, a.chr);
	strcpy(res.name, a.name);

  return(res);
}

int remove_redundant_intervals(struct exons_list *a, int num)
{
	int i = 0, j = 0;
	struct exons_list temp;
	int lo = 0, hi = 0;
	
	temp = assign_exons(a[i]);

	for ( i = 1; i < num; i++ ) 
	{
		if( ( strcmp(temp.chr, a[i].chr) == 0 ) && ( proper_overlap( temp.reg, a[i].reg ) == true ) ) {
			if( temp.reg.lower < a[i].reg.lower ) lo = temp.reg.lower;
			else lo = a[i].reg.lower;

			if( temp.reg.upper < a[i].reg.upper ) hi = a[i].reg.upper;
			else hi = temp.reg.upper;

			temp.reg = assign_I(lo, hi);
		}	
		else {
			a[j] = assign_exons(a[i-1]);
			a[j].reg = assign_I(temp.reg.lower, temp.reg.upper);
			j++;
			temp = assign_exons(a[i]);
		}
	}

	if( temp.reg.lower < a[i].reg.lower ) lo = temp.reg.lower;
	else lo = a[i].reg.lower;

	if( temp.reg.upper < a[i].reg.upper ) hi = a[i].reg.upper;
	else hi = temp.reg.upper;
	a[j] = assign_exons(a[i-1]);
	a[j].reg = assign_I(temp.reg.lower, temp.reg.upper);
	j++;
	return(j);
}

void print_exons_list(struct exons_list *a, int num)
{
	int i = 0;

	printf("## num of exons: %d\n", num);
	for( i = 0; i < num; i++ ) 
	{
		printf("%s %d %d\n", a[i].chr, a[i].reg.lower, a[i].reg.upper);	
	}
}
