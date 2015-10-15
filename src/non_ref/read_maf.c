#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"

extern char S[BIG], T[BIG];
extern int debug_mode;

char *nucs(char *x) {
  int i;

  if ((*x != 's') || (*++x != ' '))
    fatal("expecting an s-line");
  for (i = 0; i < 5; ++i) {
    while (*++x != ' ')
      ;
    while (*++x == ' ')
      ;
  }
  if (!(strchr("ACGTNBDHKMRSVWY", toupper(*x)) || strchr("-", *x)))
    fatalf("expecting nucleotides in %s", x);
  return x;
}

int cal_pid(char *s, char *t, int ncol) {
  int i, match = 0, aligned = 0;
  float pid;

  for( i = 0; i < ncol; i++ )
  {
    if( (s[i] != '-') && (toupper(s[i]) != 'N') && (toupper(t[i]) != 'N') && (t[i] != '-') )
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = (int)((((float)(match*100))/((float)aligned)) + 0.5);
  return(pid);
}

void read_maf(char *fname, int mode, struct DotList *algns, int *num_algns, int *size1, int *size2) {
	FILE *fp;
	char *status;
	int i = 0;
	int count = 0;
	int temp;
	int a_pid;
	int b1, e1, b2, e2;
	char strand[100], len1[100], len2[100];
	char *s, *t;
	int algn_type = (SELF1 - 1);
	int j = 0;
	char name1[LEN_NAME], name2[LEN_NAME];

	strcpy(name1, "");
	strcpy(name2, "");

	fp = ckopen(fname, "r");
	if (((status = fgets(S, BIG, fp)) == NULL) || strncmp(S, "##maf version", 13))
		fatalf("%s is not a maf file", fname);
	while ((status != NULL) && (S[0] == '#')) {
		if( (mode == C_MODE) && (strncmp(S, "##maf", 5) == 0) ) {
			algn_type++;
			j = 0;
		}
		status = fgets(S, BIG, fp);
	}

	while ((status != NULL) && (strstr(S, "eof") == NULL)) {
		if(S[0] == '#') {
			if(mode == C_MODE) {
				while(S[0] == '#') {
					if ((status = fgets(S, BIG, fp)) == NULL)
						fatalf("no alignments in %s", fname);
					if( strncmp(S, "##maf", 5) == 0 ) algn_type++;
				}	
				algn_type++;
				if( algn_type > PAIR ) fatal("too many alignments are combined\n");
			}
			else {
				fatalf("not supported maf format\n");
			}
			j = 0;
		}

		if (S[0] != 'a')
			fatalf("expecting an a-line in %s, saw %s",
			  fname, S);
		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
			fatalf("cannot find alignment in %s", fname);
		if ((sscanf(S, "%*s %s %d %d %*s %s", name1, &b1, &e1, len1) != 4) || (sscanf(T, "%*s %s %d %d %s %s", name2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info of 2 in %s", fname);
		}
		// aligned interval given as base-0 start and length
		e1 += b1;
		e2 += b2;

		if( strcmp(strand, "-") == 0) {
			temp = b2;
			b2 = atoi(len2) - e2;
			e2 = atoi(len2) - temp;	
		}			

		b1++;
		b2++;
//		e1++;
//		e2++;

		s = nucs(S);
		t = nucs(T);
		a_pid = cal_pid(s, t, strlen(s)-1);

		if( ((mode == D_MODE) || ((mode == C_MODE) && (algn_type <= PAIR))) && (( (algn_type != PAIR) && (b1 >= b2)) || ((algn_type != PAIR) && (abs(b1-b2) <= DEL_TH) && (abs(e1-e2) <=DEL_TH)) || ((e1-b1) < ALT_EFFEC_VALUE) || (a_pid <= PID_TH) )) {}
		else  {
			strcpy(algns[count].name1, name1);
			strcpy(algns[count].name2, name2);
			algns[count].len1 = atoi(len1);
			algns[count].len2 = atoi(len2);
			algns[count].x = assign_I(b1, e1);
			if( b2 < e2 ) algns[count].y = assign_I(b2, e2);
			else algns[count].y = assign_I(e2, b2);
			algns[count].identity = a_pid;

			if( strcmp(strand, "+") == 0 ) {
				algns[count].sign = 0;
				algns[count].init_sign = 0;
			}	
			else if( strcmp(strand, "-") == 0 ) {
				algns[count].sign = 1;
				algns[count].init_sign = 1;
			}
			else {
				algns[count].sign = DELETED;
				algns[count].init_sign = DELETED;
			}

			algns[count].fid = i; // ith alignment
			algns[count].index = count; // ith alignment
      algns[count].c_id = -1; // not chained alignment
      algns[count].m_id = -1; // not chained alignment
      algns[count].rp1_id = -1; // the inserted repeat id of the chained alignment in first seq
      algns[count].rp2_id = -1; // the inserted repeat id of the chained alignment in second seq 
      algns[count].l_id = -1;
      algns[count].lock = -1;  
      algns[count].m_x = assign_I(0,1);
      algns[count].m_y = assign_I(0,1);
      algns[count].xl_diff = 0; // the offset of the left end
      algns[count].yl_diff = 0; // the offset of the left end
      algns[count].xr_diff = 0; // the offset of the right end
      algns[count].yr_diff = 0; // the offset of the right end
      algns[count].pair_self = -1;
      algns[count].l_pid = -1;
			algns[count].sp_id = algn_type; // SELF1 for first self-alignment, SELF2 for second self-alignment and PAIR for pairwise alignment
      algns[count].xl_offset = 0; // the offset of low of x
      algns[count].yl_offset = 0; // the offset of up of x
      algns[count].xr_offset = 0; // the offset of low of y 
			if( algn_type == PAIR ) algns[count].pair_self = PAIR;
			else algns[count].pair_self = SELF;

			count++;
		}

		if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
			fatalf("bad alignment end in %s", fname);
		status = fgets(S, BIG, fp);
		i++; // ith alignment 
		j++;
	}

	*size1 = atoi(len1);
	*size2 = atoi(len2);
	*num_algns = count;
	fclose(fp);
}

float cal_pid_maf(char *s, char *t, int ncol) {
  int i, match = 0, aligned = 0;
  float pid;

  for( i = 0; i < ncol; i++ )
  {
    if( (strchr("ACGT", toupper(s[i]))) && (strchr("ACGT", toupper(t[i]))))
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = ((float)(match*100))/((float)aligned);
  return(pid);
}

float cal_pid_maf_beg(char *s, char *t, int beg, int ncol) {
  int i, match = 0, aligned = 0;
  float pid;

  for( i = beg; i < ncol; i++ )
  {
    if( (strchr("ACGT", toupper(s[i]))) && (strchr("ACGT", toupper(t[i]))))
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = ((float)(match*100))/((float)aligned);
  return(pid);
}
