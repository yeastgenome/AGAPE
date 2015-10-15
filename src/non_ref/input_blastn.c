#include "input_blastn.h"
#include "util.h"
#include "tokens.h"
#include "main.h"

int count_genes_blastn(FILE *f, int *max_hits)
{
  char buf[10000];
  int num_genes = 0;
  int num_hits = 0;
	int max_num = 0;

  while( fgets(buf, 10000, f) ) {
    if( (buf[0] == '#') && (strstr(buf, "hits found") != NULL ) ) {
      if( sscanf(buf, "%*s %d %*s", &num_hits) != 1 ) {
        fatalf("wrong blastn output in %s", buf);
      }
      else if( num_hits > 0 ) {
        num_genes++;
				if( num_hits > max_num ) {
					max_num = num_hits;
				}
      }
    }
  }
	*max_hits = max_num;
  return(num_genes);
}

bool is_name_already_in(char *name, struct sp_list *list, int num)
{
	int i = 0;
	bool res = false;

	while( (i < num) && ( res == false ) ) {
		if( strcmp(name, list[i].name) == 0 ) res = true;
		i++;
	}
	return(res);
}

int input_genes_blastn(FILE *f, struct sp_list **genes, int pid_cutoff, bool is_print)
{
  char buf[10000];
  bool is_done = false;
  char name[LEN_NAME], evalue_str[LEN_NAME], cur_name[LEN_NAME];
  float cur_pid = (float) -1;
	double evalue = (double) 1;
  int b = 0, e = 0;
  int cur_b = 0, cur_e = 0;
  int len = 0, cur_len = 0, min_len = 0;
  int num_hits = 0;
  int count = 0;
  int i = 0;
	int num_genes = 0;

  strcpy(evalue_str, "");
  strcpy(buf, "");
  strcpy(name, "");
  strcpy(cur_name, "");

	fseek(f, 0, SEEK_SET);
  while( (is_done == false) && fgets(buf, 10000, f) ) {
    while( (is_done == false) && (buf[0] == '#') ) {
      if( strstr(buf, "Query") != NULL ) {
        if( sscanf(buf, "%*s %*s %s", name) != 1 ) {
          fatalf("wrong format in %s", buf);
        }
				else {
					if( is_print == true ) printf("%s ", name);
				}
      }
      else if( strstr(buf, "hit") != NULL ) {
        if( sscanf(buf, "%*s %d %*s", &num_hits) != 1 ) {
          fatalf("wrong format in %s", buf);
        }
				else {
					if( num_hits == 0 ) {
						if( is_print == true ) printf("0\n");
					}
				}
      }
      if( !fgets(buf, 10000, f) ) is_done = true;
    }

    count = 0;
    while((is_done == false) && (buf[0] != '#')) {
      if( sscanf(buf, "%s %f %s %d %d %d %d %d %d %*s", cur_name, &cur_pid, evalue_str, &len, &b, &e, &cur_len, &cur_b, &cur_e) != 9 ) {
        fatalf("wrong format in %s", buf);
      }
			evalue = AtoF(evalue_str);

			if( len > cur_len ) min_len = cur_len;
			else min_len = len;

      cur_len = abs(cur_e - cur_b + 1);
      if( ((((float)cur_len)/((float)min_len)) > BLASTN_LEN_TH) && (evalue < EVAL_CUTOFF) && (cur_pid > pid_cutoff) )  
      {
				if( count == 0 ) {
					strcpy(genes[i][count].name, name);
					count++;
				}

				if( (count > 0 ) && (strcmp(name, cur_name) != 0) && (is_name_already_in(cur_name, genes[i], count) == false) ) {
						strcpy(genes[i][count].name, cur_name);
						count++;
        }
      }
      if( !fgets(buf, 10000, f) )is_done = true;
    }
		
		if( is_print == true ) {
			if( count > 0 ) {
				printf("%d\n", count-1);
			}
			else printf("0\n");
		}

    if( (count > 0) && (num_hits > 0) ) {
			genes[i][0].id = count;
      i++;
    }
  }

  num_genes = i;
  return(num_genes);
}
