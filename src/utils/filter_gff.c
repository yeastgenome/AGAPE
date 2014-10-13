#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_genes.h"
#include "tokens.h"

#define MIN_BASES 15
#define MIN_ORF_BASES 300

void filter_gff_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, int type);

int main(int argc, char *argv[]) {
	FILE *f;
//	int b = 0, e = 0;
	char buf[10000];
//	int sum = 0;
	struct g_list *genes1;
	struct exons_list *exons1;
	int num_genes1 = 0;
	int num_exons1 = 0;	
	int num = 0;
	int i = 0, j = 0;
	int type = 0;

	strcpy(buf, "");

	if( argc == 3 ) {
		if( strcmp(argv[2], "SGD") == 0 ) type = SGD;
		else if( strcmp(argv[2], "MAKER") == 0 ) type = MAKER;
		else if( strcmp(argv[2], "MULTI_CDS") == 0 ) type = MULTI_CDS;
		else type = STATS;
	}
	else if( argc != 2 ) {
		fatal("filter_gff gff (type)\n");
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
	
	num = input_genes_in_gff(f, genes1, exons1);	
	if( num != num_genes1 ) {
		fatalf("gene counter error in %s\n", argv[1]);
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

	filter_gff_lists(genes1, num_genes1, exons1, num_exons1, type);
	num_genes1 = rm_redun_genes(genes1, 0, num_genes1-1);
	write_in_gff(genes1, num_genes1, exons1, num_exons1);

	if( num_genes1 > 0 ) {
		free(genes1);
		free(exons1);
	}

	fclose(f);
	return EXIT_SUCCESS;
}

void filter_gff_lists(struct g_list *genes1, int num_genes1, struct exons_list *exons1, int num_exons1, int type)
{
	int i = 0, j = 0;
	struct I cur, tmp;
	int sid = 0, eid = 0;

	cur = assign_I(0, 1);
	tmp = assign_I(0, 1);

	if( type == SGD ) {
		for( i = 0; i < num_genes1; i++ ) {
    sid = genes1[i].cdsStart;
    eid = genes1[i].cdsEnd;
    if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
      cur = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
    }
    else {
      if( genes1[i].strand == '-' ) {
        cur = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
      }
      else {
        fatalf("check exons list for %s,%s:%d-%d\n", genes1[i].gname, genes1[i].sname, genes1[i].txStart, genes1[i].txEnd);
      }
    }

//			cur = assign_I(genes1[i].txStart, genes1[i].txEnd);
//			if( (width(cur) < MIN_ORF_BASES) && (strstr(genes1[i].info, "Dubious") != 0) )
			if( (genes1[i].txStart <= 0)  || (genes1[i].txEnd <= 0) ) {
				genes1[i].type = REDUN;	
			}
			else if( width(cur) < MIN_ORF_BASES )
			{
				genes1[i].type = REDUN;
			}
			else if( genes1[i].type == REDUN ) {}
			else {
				j = i+1;
				if( j < num_genes1 ) {
		      sid = genes1[j].cdsStart;
   		 		eid = genes1[j].cdsEnd;
     	 		if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
        		tmp = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
					}
      		else {
        		if( genes1[j].strand == '-' ) {
          		tmp = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
        		}
        		else {
          		fatalf("check exons list for %s,%s:%d-%d\n", genes1[j].gname, genes1[j].sname, genes1[i].txStart, genes1[j].txEnd);
        		}
      		}

//					tmp = assign_I(genes1[j].txStart, genes1[j].txEnd);
				}

				while( (j < num_genes1) && (proper_overlap(cur, tmp) == true) ) {
					if( width(tmp) < MIN_ORF_BASES )
					{
						genes1[j].type = REDUN;
					}
					else if( genes1[j].type == REDUN ) {}
					else {
						if( width(intersect(cur, tmp)) >= MIN_BASES ) {
							if( (strstr(genes1[i].info, "Verified") != 0) || (strstr(genes1[i].info, "Uncharacterized") != 0) ) {
								if(strstr(genes1[j].info, "Dubious") != 0 ) {
//									if( genes1[j].strand == genes1[i].strand ) {
										genes1[j].type = REDUN;
//									}
								}
							}
							else if( strstr(genes1[i].info, "Dubious") != 0 ) {
								if( (strstr(genes1[j].info, "Verified") != 0) || (strstr(genes1[j].info, "Uncharacterized") != 0) ) {
//									if( genes1[j].strand == genes1[i].strand ) {
										genes1[i].type = REDUN;
//									}
								}
								else if( strstr(genes1[j].info, "Dubious") != 0 ) {
									if(width(tmp) < width(cur)) {
//										if( genes1[j].strand == genes1[i].strand ) {
											genes1[j].type = REDUN;
//										}
									}
									else if(width(tmp) >= width(cur)) {
//										if( genes1[j].strand == genes1[i].strand ) {
											genes1[i].type = REDUN;
//										}
									}
								}
							}
						}
					}
	
					j++;	
					if( j < num_genes1 ) {
		   	  	sid = genes1[j].cdsStart;
   		 			eid = genes1[j].cdsEnd;
     	 			if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
       		 		tmp = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
						}
      			else {
       		 		if( genes1[j].strand == '-' ) {
       		   		tmp = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
       		 		}
       		 		else {
       		   		fatalf("check exons list for %s,%s:%d-%d\n", genes1[j].gname, genes1[j].sname, genes1[i].txStart, genes1[j].txEnd);
       		 		}
      			}
//							tmp = assign_I(genes1[j].txStart, genes1[j].txEnd);
					}
				}
			}
		}	
	}
	else if( type == MAKER ) {
		for( i = 0; i < num_genes1; i++ ) {
 	  	sid = genes1[i].cdsStart;
   		eid = genes1[i].cdsEnd;
    	if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
      	cur = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
    	}
    	else {
      	if( genes1[i].strand == '-' ) {
        	cur = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
      	}
      	else {
        	fatalf("check exons list for %s,%s:%d-%d\n", genes1[i].gname, genes1[i].sname, genes1[i].txStart, genes1[i].txEnd);
      	}
    	}
//			cur = assign_I(genes1[i].txStart, genes1[i].txEnd);
//			if( (width(cur) < MIN_ORF_BASES) && (strcmp(genes1[i].gname, "UNDEF") == 0)  )
			if( (genes1[i].type == REDUN) || (genes1[i].type == MATCH) || (genes1[i].type == PARTIAL) ) {
				genes1[i].type = REDUN;
			}
			else if( width(cur) < MIN_ORF_BASES ) 
			{
				genes1[i].type = REDUN;
			}
			else {
				j = i+1;
				if( j < num_genes1 ) {
		     	sid = genes1[j].cdsStart;
   		 		eid = genes1[j].cdsEnd;
     	 		if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
       	 		tmp = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
					}
      		else {
       	 		if( genes1[j].strand == '-' ) {
       	   		tmp = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
       	 		}
       	 		else {
       	   		fatalf("check exons list for %s,%s:%d-%d\n", genes1[j].gname, genes1[j].sname, genes1[i].txStart, genes1[j].txEnd);
       	 		}
      		}
//					tmp = assign_I(genes1[j].txStart, genes1[j].txEnd);
				}

				while( (j < num_genes1) && (proper_overlap(cur, tmp) == true) ) {
					if( (genes1[j].type == REDUN) || (genes1[j].type == MATCH) || (genes1[j].type == PARTIAL) ) {
						genes1[j].type = REDUN;
					}
					else if( width(cur) < MIN_ORF_BASES ) 
					{
						genes1[j].type = REDUN;
					}
					else {
//						if( (width(intersect(cur, tmp)) >= MIN_BASES) && (genes1[i].strand == genes1[j].strand) ) {
						if( width(intersect(cur, tmp)) >= MIN_BASES ) {
							if(width(tmp) < width(cur)) {
								genes1[j].type = REDUN;
							}
							else if(width(tmp) >= width(cur)) {
								genes1[i].type = REDUN;
							}
						}
					}
					j++;
					if( j < num_genes1 ) {
		   	  	sid = genes1[j].cdsStart;
   		 			eid = genes1[j].cdsEnd;
     	 			if( exons1[sid].reg.lower < exons1[eid].reg.upper ) {
       		 		tmp = assign_I(exons1[sid].reg.lower, exons1[eid].reg.upper);
						}
      			else {
       		 		if( genes1[j].strand == '-' ) {
       		   		tmp = assign_I(exons1[eid].reg.lower, exons1[sid].reg.upper);
       		 		}
       		 		else {
       		   		fatalf("check exons list for %s,%s:%d-%d\n", genes1[j].gname, genes1[j].sname, genes1[i].txStart, genes1[j].txEnd);
       		 		}
      			}
//						tmp = assign_I(genes1[j].txStart, genes1[j].txEnd);
					}
				}
			}
		}
	}
	else if ( type == MULTI_CDS ) {
		for( i = 0; i < num_genes1; i++ ) {
			if( genes1[i].exonCount >= 2 ) {
			}
			else {
				genes1[i].type = REDUN;
			}
		}
	}
	else {
		fatalf("Unsupported type: %d\n", type);
	}

}
