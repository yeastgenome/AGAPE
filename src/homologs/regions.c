#include "util.h"
#include "util_i.h"
#include "regions.h"

bool strict_subset(struct I reg1, struct I reg2)
{
	struct I temp;
	int len;

	if( width(reg1) > width(reg2) ) len = width(reg2);
	else len = width(reg1);

	if( len <= MIN_LEN )
	{
		temp = assign_I(reg2.lower - ERR_SM_TH, reg2.upper + ERR_SM_TH );
	}
	else if( len <= DIS_THRESHOLD )
	{
		temp = assign_I(reg2.lower - ERR_LG_TH, reg2.upper + ERR_LG_TH );
	}
	else
	{
		temp = assign_I(reg2.lower - ERR_TH, reg2.upper + ERR_TH );
	}

	return(subset(reg1, temp));
}

bool strict_almost_equal(struct I reg1, struct I reg2)
{
	if( strict_subset(reg1, reg2) == true )
	{
		if( strict_subset(reg2, reg1) == true )
		{
			return true;
		}
		else return(false);
	}
	else return(false);
}

bool almost_equal(struct I reg1, struct I reg2)
{
	int len;

	if( width(reg1) > width(reg2) ) len = width(reg2);
	else len = width(reg1);

	if( (float)((float)abs(width(reg1)-width(reg2)))/((float)len) < 0.12 )
	{
		if( proper_overlap(reg1, reg2) == true )
		{
			if( ((float)width(intersect(reg1, reg2)))/((float)len) >= 0.88 )
			{
				return true;
			}
			else return false;
		}
		else
		{
			return false;
		}
	}
	else return false;
}

bool proper_in(int r, struct I reg)
{
	if(in(r, reg)) 
	{
		if(r == reg.lower || r == reg.upper) return false;
		else return true;
	}
	else return false;
}

bool f_loose_overlap(struct I reg1, struct I reg2, int th)
{
	struct I temp1, temp2;
	int threshold;	
	int cur_th;

	cur_th = th * UNIT_THRESHOLD;
	if( (width(reg1) >= DIS_THRESHOLD) && (width(reg2) >= DIS_THRESHOLD) )
	{
		if( (3*cur_th) > D_OP_TH ) threshold = 3*cur_th;
		else threshold = D_OP_TH;

		temp1 = assign_I((reg1.lower + threshold), (reg1.upper - threshold));
		temp2 = assign_I((reg2.lower + threshold), (reg2.upper - threshold));
	}
	else if( ( width(reg1) >= (10*cur_th)) && (width(reg2) >= (10*cur_th)) )
	{
		if( (2*cur_th) > D_OP_TH ) threshold = 2*cur_th;
		else threshold = D_OP_TH;

		temp1 = assign_I((reg1.lower + threshold), (reg1.upper - threshold));
		temp2 = assign_I((reg2.lower + threshold), (reg2.upper - threshold));
	}
	else if( (width(reg1) >= (3*D_OP_TH)) && (width(reg2) >= (3*D_OP_TH)) )
	{
		threshold = D_OP_TH;

		temp1 = assign_I((reg1.lower + threshold), (reg1.upper - threshold));
		temp2 = assign_I((reg2.lower + threshold), (reg2.upper - threshold));
	}
	else if( (width(reg1) >= (6*UNIT_THRESHOLD)) && (width(reg2) >= (6*UNIT_THRESHOLD)))
	{
		temp1 = assign_I((reg1.lower + UNIT_THRESHOLD), (reg1.upper - UNIT_THRESHOLD));
		temp2 = assign_I((reg2.lower + UNIT_THRESHOLD), (reg2.upper - UNIT_THRESHOLD));
	}
	else if( (width(reg1) >= (5*SMALLEST_TH)) && (width(reg2) >= (5*SMALLEST_TH)))
	{
		temp1 = assign_I((reg1.lower + SMALLEST_TH), (reg1.upper - SMALLEST_TH));
		temp2 = assign_I((reg2.lower + SMALLEST_TH), (reg2.upper - SMALLEST_TH));
	}
	else 
	{
		temp1 = assign_I(reg1.lower, reg1.upper);
		temp2 = assign_I(reg2.lower, reg2.upper);
	}

	return(proper_overlap(temp1, temp2)); 
}

bool loose_overlap(struct I reg1, struct I reg2)
{
	struct I temp1, temp2;

	if( (width(reg1) >= DIS_THRESHOLD) && (width(reg2) >= DIS_THRESHOLD) )
	{
		temp1 = assign_I((reg1.lower + 3*THRESHOLD), (reg1.upper - 3*THRESHOLD));
		temp2 = assign_I((reg2.lower + 3*THRESHOLD), (reg2.upper - 3*THRESHOLD));
	}
	else if( ( width(reg1) >= (10*THRESHOLD) ) && (width(reg2) >= (10*THRESHOLD)) )
	{
		temp1 = assign_I((reg1.lower + 2*THRESHOLD), (reg1.upper - 2*THRESHOLD));
		temp2 = assign_I((reg2.lower + 2*THRESHOLD), (reg2.upper - 2*THRESHOLD));
	}
	else if( (width(reg1) >= (10*SMALLEST_TH)) && (width(reg2) >= (10*SMALLEST_TH)))
	{
		temp1 = assign_I((reg1.lower + SMALLEST_TH), (reg1.upper - SMALLEST_TH));
		temp2 = assign_I((reg2.lower + SMALLEST_TH), (reg2.upper - SMALLEST_TH));
	}
	else 
	{
		temp1 = assign_I(reg1.lower, reg1.upper);
		temp2 = assign_I(reg2.lower, reg2.upper);
	}

	return(overlap(temp1, temp2)); 
}

bool tight_overlap(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg2) > (4*T_OP_TH) )
	{
		temp = assign_I((reg2.lower + T_OP_TH), (reg2.upper - T_OP_TH));
		return(overlap(reg1, temp)); 
	}
	else 
	{
		return(overlap(reg1, reg2));
	}
}

bool proper_overlap(struct I reg1, struct I reg2)
{
	if(overlap(reg1, reg2)) 
	{
		if(width(intersect(reg1, reg2)) == 0) return false;
		else return true;
	}
	else return false;
}

bool strict_overlap(struct I reg1, struct I reg2, int threshold)
{
	struct I temp;

	temp = assign_I((reg2.lower - threshold), (reg2.upper + threshold));

	return(overlap(reg1, temp)); 
}

bool left_loose_overlap(struct I reg1, struct I reg2, int th)
{
	struct I temp;
	int threshold;	

	threshold = th * UNIT_THRESHOLD;
	if( width(reg2) >= DIS_THRESHOLD )
	{
		temp = assign_I((reg2.lower + (3*threshold)), (reg2.upper - (3*threshold)));
	}
	else if( width(reg2) >= (10*threshold) )
	{
		temp = assign_I((reg2.lower + (2*threshold)), (reg2.upper - (2*threshold)));
	}
	else if( width(reg2) >= (4*UNIT_THRESHOLD))
	{
		temp = assign_I((reg2.lower + UNIT_THRESHOLD), (reg2.upper - UNIT_THRESHOLD));
	}
	else if( width(reg2) >= (10*SMALLEST_TH) )
	{
		temp = assign_I((reg2.lower + SMALLEST_TH), (reg2.upper - SMALLEST_TH));
	}
	else 
	{
		temp = assign_I(reg2.lower, reg2.upper);
	}

	return(in(reg1.lower, temp)); 
}

bool right_loose_overlap(struct I reg1, struct I reg2, int th)
{
	struct I temp;
	int threshold;	

	threshold = th * UNIT_THRESHOLD;
	if( width(reg2) >= DIS_THRESHOLD )
	{
		temp = assign_I((reg2.lower + (3*threshold)), (reg2.upper - (3*threshold)));
	}
	else if( width(reg2) >= (10*threshold) )
	{
		temp = assign_I((reg2.lower + (2*threshold)), (reg2.upper - (2*threshold)));
	}
	else if( width(reg2) >= (4*UNIT_THRESHOLD))
	{
		temp = assign_I((reg2.lower + UNIT_THRESHOLD), (reg2.upper - UNIT_THRESHOLD));
	}
	else if( width(reg2) >= (10*SMALLEST_TH) )
	{
		temp = assign_I((reg2.lower + SMALLEST_TH), (reg2.upper - SMALLEST_TH));
	}
	else 
	{
		temp = assign_I(reg2.lower, reg2.upper);
	}

	return(in(reg1.upper, temp)); 
}

bool fully_subset(struct I reg1, struct I reg2)
{
	struct I temp;
	
	if( width(reg2) > 8*THRESHOLD )
	{
		temp = assign_I((reg2.lower) + 4*THRESHOLD, (reg2.upper) - 4*THRESHOLD);
	}
	else if( width(reg2) > 4*THRESHOLD )
	{
		temp = assign_I((reg2.lower) + 2*THRESHOLD, (reg2.upper) - 2*THRESHOLD);
	}
	else temp = assign_I(reg2.lower, reg2.upper);
	
	return(subset(reg1, temp));
}

bool almost_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	temp = assign_I((reg2.lower) - THRESHOLD, (reg2.upper) + THRESHOLD);
	if( (width(reg1) > 5*THRESHOLD) && (width(reg2) > 5*THRESHOLD))
	{
	  if(subset(reg1, temp))
		{
			return(true);
		}
		else if( proper_overlap(reg1, reg2) == true )
		{
			if(width(intersect(reg1, reg2)) > ((0.8)*width(reg1)))
				return(true);
			else return(false);
		}
		else return(false);
	}
	else
	{
		return(subset(reg1, temp));
	}
}

bool f_loosen_subset(struct I reg1, struct I reg2, int th)
{
	struct I temp;
	int threshold;

	threshold = th * UNIT_THRESHOLD;

	if( width(reg1) > (DIS_THRESHOLD) ) 
	{
		temp = assign_I((reg2.lower) - (5*threshold), (reg2.upper) + (5*threshold));
	}
	else if( width(reg1) > (DIS_THRESHOLD/2) ) 
	{
		temp = assign_I((reg2.lower) - (3*threshold), (reg2.upper) + (3*threshold));
	}
	else
	{
		temp = assign_I((reg2.lower) - (2*threshold), (reg2.upper) + (2*threshold));
	}

	return(subset(reg1, temp));
}

bool loosen_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg1) > (DIS_THRESHOLD) ) 
	{
		temp = assign_I((reg2.lower) - (5*THRESHOLD), (reg2.upper) + (5*THRESHOLD));
	}
	else
	{
		temp = assign_I((reg2.lower) - (2*THRESHOLD), (reg2.upper) + (2*THRESHOLD));
	}
	return(subset(reg1, temp));
}

bool too_loosen_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg1) > LOOSEN_T ) 
	{
		temp = assign_I((reg1.lower) + (LOOSEN_T/2), (reg1.upper) - (LOOSEN_T/2));
	}
	else
	{
		if( width(reg1) > (LOOSEN_T/2) )
		{	
			temp = assign_I((reg1.lower) + (LOOSEN_T/4), (reg1.upper) - (LOOSEN_T/4));
		}
		else
		{
			fatal("error: less than DIS_THRESHOLD\n");
		}
	}
	return(subset(temp, reg2));
}

void init_array(int *array, int num)
{
	int i;

	for( i = 0 ; i < num; i++ ) array[i] = 0;
}

void overwrite_dots(int *num, struct DotList *dots)
{
  int i;
	int j = 0;

	for(i = 0 ; i < *num; i++)
	{
		if( dots[i].l_id != -1 )
		{
			dots[dots[i].l_id].sign = dots[i].sign;
			dots[i].x = dots[i].m_x;
			dots[i].y = dots[i].m_y;
			dots[i].l_id = -1;
			dots[i].m_x = assign_I(0,1);
			dots[i].m_y = assign_I(0,1);
			dots[i].identity = dots[i].m_pid;
			dots[i].m_pid = 0;
		}
	}

  for(i = 0 ; i < *num; i++)
  {
    if( (dots[i].sign == 0) || (dots[i].sign == 1) || (dots[i].sign == ORTHO) || (dots[i].sign == ORTHO_COMP) )
    {
      dots[j].x = assign_I(dots[i].x.lower, dots[i].x.upper);
      dots[j].y = assign_I(dots[i].y.lower, dots[i].y.upper);
	    dots[j].sign = dots[i].sign;
	    dots[j].identity = dots[i].identity;
			dots[j].l_id = dots[i].l_id;
			dots[j].lock = -1;
			dots[j].c_id = dots[i].c_id;
			dots[j].m_id = dots[i].m_id;
			dots[j].fid = dots[i].fid;
			dots[j].index = dots[i].index;
			dots[j].l_pid = dots[i].l_pid;
			dots[j].xl_diff = dots[i].xl_diff;
			dots[j].xr_diff = dots[i].xr_diff;
			dots[j].yl_diff = dots[i].yl_diff;
			dots[j].yr_diff = dots[i].yr_diff;
			dots[j].m_x = assign_I(0, 1);
			dots[j].m_y = assign_I(0, 1);
			dots[j].m_pid = 0;
			dots[j].pair_self = dots[i].pair_self;
			dots[j].sp_id = dots[i].sp_id;
			dots[j].rp1_id = dots[i].rp1_id;
			dots[j].rp2_id = dots[i].rp2_id;
      j++;
    }
		else 
		{
		}
  }
  *num = j; // Modify the number of lines
}

int compute_distance(struct I x1, struct I y1, struct I x2, struct I y2, int sign)
{   
  int c, d;
  int min;
  int x1_m, y1_m, x2_m, y2_m;

  x1_m = (x1.lower + x1.upper)/2; 
  y1_m = (y1.lower + y1.upper)/2;
  x2_m = (x2.lower + x2.upper)/2;
  y2_m = (y2.lower + y2.upper)/2;
    
  if( sign == 0 )
  {
    c = y2.lower - x2.lower;
    d = abs(x1.lower - y1.lower + c);
    min = d;
    d = abs(x1.upper - y1.upper + c);
    if( min > d ) min = d;
    d = abs(x1_m - y1_m + c);
    if( min > d ) min = d;
  
    c = y2.upper - x2.upper;
    d = abs(x1.lower - y1.lower + c);
    if( min > d ) min = d;
    d = abs(x1.upper - y1.upper + c);
    if( min > d ) min = d;
    d = abs(x1_m - y1_m + c);
    if( min > d ) min = d;
  
    c = y2_m - x2_m;
    d = abs(x1.lower - y1.lower + c);
    if( min > d ) min = d;
    d = abs(x1.upper - y1.upper + c);
    if( min > d ) min = d;
    d = abs(x1_m - y1_m + c);
    if( min > d ) min = d;
  }
  else
  {
    c = (-1)*(x2.lower + y2.upper);
    d = abs(x1.lower + y1.upper + c);
    min = d;
    d = abs(x1.upper + y1.lower + c);
    if( min > d ) min = d;
    d = abs(x1_m + y1_m + c);
    if( min > d ) min = d;

    c = (-1)*(x2.upper + y2.lower);
    d = abs(x1.lower + y1.upper + c);
    if( min > d ) min = d;
    d = abs(x1.upper + y1.lower + c);
    if( min > d ) min = d;
    d = abs(x1_m + y1_m + c);
    if( min > d ) min = d;

    c = (-1)*(x2_m + y2_m);
    d = abs(x1.lower - y1.lower + c);
    if( min > d ) min = d;
    d = abs(x1.upper - y1.upper + c);
    if( min > d ) min = d;
    d = abs(x1_m - y1_m + c);
    if( min > d ) min = d;
  }

  return min;
}
