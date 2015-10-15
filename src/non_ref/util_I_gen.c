#include "main.h"
#include "util_I_gen.h"
#include "util_i.h"
#include "util_gen.h"
#include "util.h"

int sort_merge_intervals_and_pid(struct I *regs, int num_regs, struct I *new_regs, int *pid) // closed interval
{
  int j = 0;
  int i = 0;
  int lo = 0, hi = 0;
	int pid_sum = 0;
	int total_bases = 0;

  if( num_regs > 0 ) {
    quick_sort_inc_int(regs, 0, num_regs-1);
    j =  0;
    new_regs[j] = assign_I(regs[0].lower, regs[0].upper);
		pid_sum = width(regs[0]) * pid[0];
		total_bases = width(regs[0]);
    for( i = 1; i < num_regs; i++ ) {
      while( (i < num_regs) && (overlap(new_regs[j], regs[i]) == true ) ) {
        if( new_regs[j].lower < regs[i].lower ) lo = new_regs[j].lower;
        else lo = regs[i].lower;

        if( new_regs[j].upper < regs[i].upper ) hi = regs[i].upper;
        else hi = new_regs[j].upper;

				pid_sum = pid_sum + (width(regs[i]) * pid[i]);
				total_bases = total_bases + width(regs[i]);
        new_regs[j] = assign_I(lo, hi);
        i++;
      }

      if( i >= num_regs ) {
			}
      else {
				pid[j] = (int)(((float)pid_sum / (float)total_bases) + 0.5);
        j++;
        if( j > num_regs ) {
          fatalf("overflow in new_regs[]: %d\n", j);
        }
        new_regs[j] = assign_I(regs[i].lower, regs[i].upper);
				pid_sum = width(regs[i]) * pid[i];
				total_bases = width(regs[i]);
      }
    }

    j++;
  }

	for( i = j; i < num_regs; i++ ) pid[i] = 0;
	
  return(j);
}

int sort_merge_intervals(struct I *regs, int num_regs, struct I *new_regs) // closed interval
{
  int j = 0;
  int i = 0;
  int lo = 0, hi = 0;

  if( num_regs > 0 ) {
    quick_sort_inc_int(regs, 0, num_regs-1);
    j =  0;
    new_regs[j] = assign_I(regs[0].lower, regs[0].upper);
    for( i = 1; i < num_regs; i++ ) {
      while( (i < num_regs) && (overlap(new_regs[j], regs[i]) == true ) ) {
        if( new_regs[j].lower < regs[i].lower ) lo = new_regs[j].lower;
        else lo = regs[i].lower;

        if( new_regs[j].upper < regs[i].upper ) hi = regs[i].upper;
        else hi = new_regs[j].upper;

        new_regs[j] = assign_I(lo, hi);
        i++;
      }

      if( i >= num_regs ) {}
      else {
        j++;
        if( j > num_regs ) {
          fatalf("overflow in new_regs[]: %d\n", j);
        }
        new_regs[j] = assign_I(regs[i].lower, regs[i].upper);
      }
    }

    j++;
  }
  return(j);
}

void initialize_I_list(struct I *list, int count)
{
  int i = 0;

  for( i = 0; i < count; i++ ) {
    list[i] = assign_I(0, 1);
  }
}

void initialize_orf_I_list(struct orf_I *list, int count)
{
  int i = 0;

  for( i = 0; i < count; i++ ) {
    list[i].region = assign_I(0, 1);
		strcpy(list[i].name, "");
		strcpy(list[i].strain_name, "");
  }
}

void initialize_scf_I_list(struct scf_I *list, int count)
{
  int i = 0;

  for( i = 0; i < count; i++ ) {
    list[i].region = assign_I(0, 1);
		strcpy(list[i].name, "");
  }
}

int input_orf_I_list(FILE *f, struct orf_I *match_regions, int count)
{
	int i = 0;
	char buf[MAX_NAME];
	char strain_name[MAX_NAME];
	char scf_name[MAX_NAME];
	char temp1[MAX_NAME], temp2[MAX_NAME];
	int b = 0, e = 0;
	char num1[100], num2[100];
	int num_match_regions = 0;

	strcpy(buf, "");
	strcpy(strain_name, "");
	strcpy(scf_name, "");
	strcpy(temp1, "");
	strcpy(temp2, "");
	strcpy(num1, "");
	strcpy(num2, "");

	fseek(f, 0, SEEK_SET);
	while( fgets(buf, MAX_NAME, f) ) {
		if( sscanf(buf, "%s %s", strain_name, scf_name) != 2 ) {
			fatalf("wrong line : %s", buf);
		}
		else {
			if( sscanf(strain_name, "%[^.].%*s", temp1) != 1 ) {
				fatalf("wrong ORF name : %s", strain_name);
			}
			else strcpy(strain_name, temp1);
			
			if( sscanf(scf_name, "%[^:]:%s", temp1, temp2) != 2 ) {
				fatalf("wrong scf interval : %s", scf_name);
			}
			else {
				strcpy(scf_name, temp1);
			}

			if( sscanf(temp2, "%[^-]-%s", num1, num2) != 2 ) {
				fatalf("wrong interval : %s", temp2);
			}
			else {
				b = atoi(num1);
				e = atoi(num2);
			}

			if( i >= count ) {
				fatalf("counting error while reading list: %d\n", count);
			}
			match_regions[i].region = assign_I(b, e);
			strcpy(match_regions[i].name, scf_name);
			strcpy(match_regions[i].strain_name, strain_name);
			i++;
		}
	}
	num_match_regions = i;
	return(num_match_regions);
}

int input_scf_I_list(FILE *f, struct scf_I *match_regions, int count)
{
	int i = 0;
	char buf[MAX_NAME];
	char cur_name[MAX_NAME];
	int b = 0, e = 0;
	int num_match_regions = 0;

	strcpy(buf, "");
	strcpy(cur_name, "");
	fseek(f, 0, SEEK_SET);
	while( fgets(buf, MAX_NAME, f) ) {
		if( sscanf(buf, "%s %d %d", cur_name, &b, &e) != 3 ) {
			fatalf("wrong line : %s", buf);
		}
		else {
			if( i >= count ) {
				fatalf("counting error while reading list: %d\n", count);
			}
			match_regions[i].region = assign_I(b, e);
			strcpy(match_regions[i].name, cur_name);
			i++;
		}
	}
	num_match_regions = i;
	return(num_match_regions);
}

int find_match_regions(struct I *cur_regions, int num, struct scf_I *match_regions, int num_match_regions, char *scaf_name)
{
  int i = 0, j = 0;

  for( i = 0; i < num_match_regions; i++ )
  {
    if( strcmp(match_regions[i].name, scaf_name) == 0 ) {
      cur_regions[num+j] = assign_I(match_regions[i].region.lower, match_regions[i].region.upper);
      j++;
    }
  }
  return(j);
}
