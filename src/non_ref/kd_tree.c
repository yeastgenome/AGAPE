#include "main.h"
#include "util_gen.h"
#include "util.h"
#include "kd_tree.h"

struct kdnode *build_kd(struct perm_pt *p_pts, int l, int u)
{
	int cutoff = 20;
	int m;
	struct kdnode *p;

	p = (struct kdnode *) ckalloc(sizeof(struct kdnode));

	if( u-l+1 <= cutoff )
	{
		p->bucket = 1;
		p->lopt = l;
		p->hipt = u;
		p->loson = NULL;
		p->hison = NULL;
	}
	else
	{
		p->bucket = 0;
		p->cutdim = findmaxspread(p_pts, l, u);
		if( p->cutdim == 1 ) quick_sort_plist_x(p_pts, l, u);
		else if( p->cutdim == 2 ) quick_sort_plist_y(p_pts, l, u);
		m = (l + u) / 2;
		m = adjust_m(p_pts, l, u, m, p->cutdim);
		p->cutval = px(p_pts, m, p->cutdim);

		if( p->cutdim == 1 ) 
		{
			p->loson = build_kd(p_pts, l, m-1);
			p->hison = build_kd(p_pts, m+1, u);
		}
		else
		{
			p->loson = build_kd(p_pts, l, m-1);
			p->hison = build_kd(p_pts, m+1, u);
		}
	}

	return p;
}

void free_kd(struct kdnode *t)
{
	if(t == NULL) return;
	else
	{
		if( t->loson != NULL) free_kd(t->loson);
		if( t->hison != NULL) free_kd(t->hison);
		free(t);
	}
}

int adjust_m(struct perm_pt *p_pts, int l, int u, int m, int cutdim)
{
	int i = m;
	
	if( cutdim == 1 )
	{
		while(((i+1) <= u) && ( p_pts[i+1].x_pt == p_pts[m].x_pt )) i++;

		if( i <= u ) return(i);
		else return(u);
	}	
	else
	{
		while(((i-1) >= l) && ( p_pts[i-1].x_pt == p_pts[m].x_pt )) i--;

		if( i >= l ) return(i);
		else return(l);
	}
}

int px(struct perm_pt *p_pts, int m, int cutdim)
{
	if( cutdim == 1 ) return(p_pts[m].x_pt);
	else return(p_pts[m].y_pt);
}

void assign_perm(struct perm_pt *p_pts, int num, struct DotList *points, int side)
{
	int i;

	for( i = 0; i < num; i++ )
	{
	  p_pts[i].id = i;
	  p_pts[i].sign = points[i].sign;
	  if( p_pts[i].sign == 0 )
	  {
			if( side == LEFT )
			{
				p_pts[i].x_pt = points[i].x.upper;
				p_pts[i].y_pt = points[i].y.upper;
			}
			else if( side == RIGHT )
			{
				p_pts[i].x_pt = points[i].x.lower;
				p_pts[i].y_pt = points[i].y.lower;
			}
	  }
	  else if( p_pts[i].sign == 1 )
	  {
			if( side == LEFT )
			{
				p_pts[i].x_pt = points[i].x.upper;
				p_pts[i].y_pt = points[i].y.lower;
			}
			else if( side == RIGHT )
			{
				p_pts[i].x_pt = points[i].x.lower;
				p_pts[i].y_pt = points[i].y.upper;
			}
	  }
	}
}

int findmaxspread(struct perm_pt *p_pts, int l, int u)
{
	int i;
	int x_max = 0, x_min = Max_base;
	int y_max = 0, y_min = Max_base;
	
	for( i = l; i <= u; i++ )
	{
		if( p_pts[i].x_pt < x_min ) x_min = p_pts[i].x_pt;
		if( p_pts[i].x_pt > x_max ) x_max = p_pts[i].x_pt;
		if( p_pts[i].y_pt < y_min ) y_min = p_pts[i].y_pt;
		if( p_pts[i].y_pt > y_max ) y_max = p_pts[i].y_pt;
	}

	if( (x_max - x_min) >= (y_max - y_min) ) return(1);
	else return(2);
}
