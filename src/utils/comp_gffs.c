#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"
#include "tokens.h"

void comp_gff_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, struct g_list *genes2, int num_genes2, struct exons_list *exons2, int num_exons2, int type);

int main(int argc, char *argv[]) {
	FILE *f, *g;
//	int b = 0, e = 0;
	char buf[10000];
//	int sum = 0;
	struct g_list *genes1, *genes2;
	struct exons_list *exons1, *exons2;
	int num_genes1 = 0, num_genes2 = 0;
	int num_exons1 = 0, num_exons2 = 0;	
	int num = 0;
	int i = 0, j = 0;
	int type = 0;

	strcpy(buf, "");

	if( argc == 4 ) {
		if( strcmp(argv[3], "SAME_STOP") == 0 ) type = SAME_STOP;
		else if( strcmp(argv[3], "OVERLAP") == 0 ) type = OVERLAP;
		else if( strcmp(argv[3], "MISSED") == 0 ) type = MISSED;
		else if( strcmp(argv[3], "NOVEL") == 0 ) type = NOVEL;
		else if( strcmp(argv[3], "IDENTICAL") == 0 ) type = IDENTICAL;
		else if( strcmp(argv[3], "CORRECT") == 0 ) type = CORRECT;
		else if( strcmp(argv[3], "COMPLETE_MISS") == 0 ) type = COMPLETE_MISS;
		else type = STATS;
	}
	else if( argc != 3 ) {
		fatal("merge_gff codex gff2 (type)\n");
	}
	

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		num_genes1 = count_genes_in_gff(f, &num_exons1);
		if( num_genes1 > 0 ) {
			genes1 = (struct g_list *) ckalloc(num_genes1 * sizeof(struct g_list));	
			if( num_exons1 < num_genes1 ) num_exons1 = num_genes1;
			exons1 = (struct exons_list *) ckalloc(num_exons1 * sizeof(struct exons_list));	
			initialize_genes(genes1, num_genes1);
			initialize_exons(exons1, num_exons1);
		}
	}
	fseek(f, 0, SEEK_SET);
	
	if((g = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[2]);
	}
	else {
		num_genes2 = count_genes_in_gff(g, &num_exons2);
		if( num_genes2 > 0 ) {
			genes2 = (struct g_list *) ckalloc(num_genes2 * sizeof(struct g_list));	
			if( num_exons2 < num_genes2 ) num_exons2 = num_genes2;
			exons2 = (struct exons_list *) ckalloc(num_exons2 * sizeof(struct exons_list));	
			initialize_genes(genes2, num_genes2);
			initialize_exons(exons2, num_exons2);
		}
	}
	fseek(g, 0, SEEK_SET);

	num = input_genes_in_gff(f, genes1, exons1);	
	if( num != num_genes1 ) {
		fatalf("gene counter error in %s\n", argv[1]);
	}

	num = input_genes_in_gff(g, genes2, exons2);	
	if( num != num_genes2 ) {
		fatalf("gene counter error in %s\n", argv[2]);
	}

	if( num_genes1 > 0 ) {
		quick_sort_inc_genes(genes1, 0, num_genes1-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes1 ) {
		j = 0;
    while( ((i+j) < num_genes1) && (genes1[i].txStart == genes1[i+j].txStart )) j++;
    quick_sort_dec_genes(genes1, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	if( num_genes2 > 0 ) {
		quick_sort_inc_genes(genes2, 0, num_genes2-1, POS_BASE);
	}
	i = 0;
	while( i < num_genes2 ) {
		j = 0;
    while( ((i+j) < num_genes2) && (genes2[i].txStart == genes2[i+j].txStart )) j++;
    quick_sort_dec_genes(genes2, i, i+j-1, LEN_BASE);
    i = i+j;
	}

	comp_gff_lists(genes1, num_genes1, exons1, num_exons1, genes2, num_genes2, exons2, num_exons2, type);

	if( num_genes1 > 0 ) {
		free(genes1);
		free(exons1);
	}

	if( num_genes2 > 0 ) {
		free(genes2);
		free(exons2);
	}

	fclose(f);
	fclose(g);
	return EXIT_SUCCESS;
}

void comp_gff_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, struct g_list *genes2, int num_genes2, struct exons_list *exons2, int num_exons2, int type)
{
	int i = 0, j = 0;
	struct I cur, tmp;
	struct I cur_CDS, tmp_CDS; 
	int sid = 0, eid = 0;
	bool is_same = false;
	bool same_stop = false;
	bool is_overlap = false;
	int num_equal = 0;
	int num_same_stop = 0;
	int num_overlap = 0;
	bool missed[num_genes2];
	bool novel[num_genes1];
	int num_missed = 0;
	int num_novel = 0;
  char gname[2][LEN_NAME];
  char seq_name[2][LEN_NAME];
  char chr[2][LEN_NAME];
  int b[2], e[2];
  float pid[2];
  char line[2][LEN_LONG_NAME];
	char SGD_gname[LEN_NAME];
	int num_complete_miss = 0;
	int complete_miss[num_genes2];
	int max_width = 0;

  for( i = 0; i < 2; i++ ) {
	  strcpy(gname[i], "");
    strcpy(seq_name[i], "");
    strcpy(chr[i], "");
    b[i] = 0;
    e[i] = 0;
    pid[i] = (float) 0;
  }
	strcpy(SGD_gname, "");

	cur = assign_I(0, 1);
	tmp = assign_I(0, 1);
	cur_CDS = assign_I(0, 1);
	tmp_CDS = assign_I(0, 1);

	for( i = 0; i < num_genes1; i++ ) {
		novel[i] = false;
	}

	for( i = 0; i < num_genes2; i++ ) {
		missed[i] = true;
		complete_miss[i] = -1;
	}

/*
	for( i = 0; i < num_genes2; i++ ) {
		cur = assign_I(genes2[i].txStart, genes2[i].txEnd);
		j = 0;
		is_overlap = false;
		while( (j < num_genes1) && (is_overlap == false) ) {
			tmp = assign_I(genes1[j].txStart, genes1[j].txEnd);
			if (proper_overlap(cur, tmp) == true) {
				is_overlap = true;
			}
			j++;
		}

		if( is_overlap == true ) {
			missed[i] = false;
		}
		else {
			j = 0;	
		}
	}
*/

	for( i = 0; i < num_genes1; i++ ) {
    if( strstr(genes1[i].gname, ";") != NULL ) {
      if( sscanf(genes1[i].gname, "%[^;];%s", line[0], line[1]) != 2 )
      {
        fatalf("wrong format in the last column: %s\n", genes1[i].gname);
      }
    }
    else {
      strcpy(line[0], genes1[i].gname);
      strcpy(line[1], genes1[i].gname);
    }

		if( strcmp(line[0], "UNDEF") == 0 ) {
			strcpy(gname[0], "UNDEF");
		}
		else {
			parse_last_column_gff(line[0], gname[0], seq_name[0], chr[0], &b[0], &e[0], &pid[0]);
		}

		if( strcmp(line[1], "UNDEF") == 0 ) {
			strcpy(gname[1], "UNDEF");
		}
		else {
    	parse_last_column_gff(line[1], gname[1], seq_name[1], chr[1], &b[1], &e[1], &pid[1]);
		}
	
		cur = assign_I(genes1[i].txStart, genes1[i].txEnd);
		sid = genes1[i].cdsStart;
		eid = genes1[i].cdsEnd;
		if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
			cur_CDS = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
		}
		else {
			if( genes1[i].strand == '-' ) {
				cur_CDS = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
			}
			else {
				fatalf("check exons list for %s,%s:%d-%d\n", genes1[i].gname, genes1[i].sname, genes1[i].txStart, genes1[i].txEnd);
			}
		}

		max_width = 0;
		for( j = 0; j < num_genes2; j++ ) {
			tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			sid = genes2[j].cdsStart;
			eid = genes2[j].cdsEnd;
			if( exons2[sid].reg.lower < exons2[eid].reg.upper ) {
				tmp_CDS = assign_I(exons2[sid].reg.lower, exons2[eid].reg.upper);
			}
			else {
				if( genes2[j].strand == '-' ) {
					tmp_CDS = assign_I(exons2[eid].reg.lower, exons2[sid].reg.upper);
				}
				else {
					fatalf("check exons list for %s,%s:%d-%d\n", genes2[j].gname, genes2[j].sname, genes2[i].txStart, genes2[j].txEnd);
				}
			}

//			if( (proper_overlap(cur, tmp) == true)  && ((width(intersect(cur, tmp)) >= ((width(cur))/2)) || (width(intersect(cur, tmp)) >= (width(tmp)/2) ) ) && ((strstr(genes2[j].gname, gname[0]) != 0 ) || (strstr(genes2[j].gname, gname[1]) != 0)) ) {
			if( (equal(cur, tmp) == true) || (equal(cur_CDS, tmp_CDS) == true) ) {
				missed[j] = false;
				if( type == CORRECT ) {
					print_item_gff(genes1[i], exons1);	
					print_item_gff(genes2[j], exons2);	
				}
			}
			
			if( (proper_overlap(cur, tmp) == true) && ((width(intersect(cur, tmp)) >= ((width(cur))/2)) || (width(intersect(cur, tmp)) >= (width(tmp)/2) ) ) ) 
			{
				if( max_width < width(intersect(cur, tmp)) ) {
					complete_miss[j] = i;
				}
			}
		}

		j = 0;
		is_same = false;
		while( (j < num_genes2) && (is_same == false) ) {
			tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			sid = genes2[j].cdsStart;
			eid = genes2[j].cdsEnd;
			if( exons2[sid].reg.lower < exons2[eid].reg.upper ) {
				tmp_CDS = assign_I(exons2[sid].reg.lower, exons2[eid].reg.upper);
			}
			else {
				if( genes2[j].strand == '-' ) {
					tmp_CDS = assign_I(exons2[eid].reg.lower, exons2[sid].reg.upper);
				}
				else {
					fatalf("check exons list for %s,%s:%d-%d\n", genes2[j].gname, genes2[j].sname, genes2[i].txStart, genes2[j].txEnd);
				}
			}

			if( (equal(cur, tmp) == true) || (equal(cur_CDS, tmp_CDS) == true) ) {
				if( type == IDENTICAL ) {
					print_item_gff(genes1[i], exons1);
				}
				is_same = true;
				num_equal++;
			}
			j++;
		}

		j = 0;
		same_stop = false;
		while( (j < num_genes2) && (is_same == false) && (same_stop == false) ) {
			tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			if( proper_overlap(cur, tmp) == true ) {
				if( genes1[i].strand == genes2[j].strand ) {
					if( genes1[i].strand == '+') {
						if( genes1[i].txEnd == genes2[j].txEnd ) {
							if( type == SAME_STOP ) {
								print_item_gff(genes1[i], exons1);
							}
							same_stop = true;
							num_same_stop++;
						}
					}
					else {
						if( genes1[i].txStart == genes2[j].txStart ) {
							if( type == SAME_STOP ) {
								print_item_gff(genes1[i], exons1);
							}
							same_stop = true;
							num_same_stop++;
						}
					}
				}
			}
			j++;
		}

		j = 0;
		is_overlap = false;
		while( (j < num_genes2) && (is_same == false) && (same_stop == false) && (is_overlap == false)) {
			tmp = assign_I(genes2[j].txStart, genes2[j].txEnd);
			if( proper_overlap(cur, tmp) == true ) {
				if( genes1[i].strand == genes2[j].strand ) {
					if( type == OVERLAP ) {
						print_item_gff(genes1[i], exons1);
					}
				}
				else {
					if( type == OVERLAP ) {
						print_item_gff(genes1[i], exons1);
					}
				}
				is_overlap = true;
				num_overlap++;
			}
			j++;
		}

		if( (is_same == false) && (same_stop == false) && (is_overlap == false) ) {
			if( type == NOVEL ) {
				print_item_gff(genes1[i], exons1);
			}
			novel[i] = true;
			num_novel++;
		}
	}

	for( i = 0; i < num_genes2; i++ ) {
		if( missed[i] == true ) {
			num_missed++;
			if( type == MISSED ) {
				print_item_gff(genes2[i], exons2);
			}
			else if( type == COMPLETE_MISS ) {
				if( complete_miss[i] == -1 ) {
					print_item_gff(genes2[i], exons2);
//					print_item_gff(genes1[complete_miss[i]], exons1);
				}
			}
		}

		if( complete_miss[i] == -1 ) {
			num_complete_miss++;
		}
	}

	if( type == STATS ) {
//		printf("%d %d %d %d %d %d %d %d\n", num_genes1, num_equal, num_same_stop, num_overlap, num_novel, num_genes2, num_missed, num_complete_miss);
		printf("%d %d %d %d %d %d %d\n", num_genes1, num_equal, num_same_stop, num_overlap, num_novel, num_genes2, num_missed);
	}
}
