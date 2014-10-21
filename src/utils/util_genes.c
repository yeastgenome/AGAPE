#include "main.h"
#include "util_genes.h"
#include "util.h"
#include "util_i.h"
#include "tokens.h"

#define OV_TH 15
#define CUT_RATIO 0.9
#define EVAL_CUTOFF 1e-06
#define MIN_EVAL_CUTOFF 1e-100
#define PID_CUTOFF 80
#define MIN_NU_LEN 300
#define MIN_LEN_RATIO 0.05
#define HIGH_PID 95

extern int debug_mode;

void initialize_scaffolds(struct scaffold *a, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) {
		a[i].b = 0;
		a[i].e = 0;
		a[i].len = 0;
		strcpy(a[i].name, "");
	}
}

void initialize_genes(struct g_list *a, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) {
		a[i].gid = -1;
		a[i].sid = -1;
		a[i].strand = '+';
		a[i].txStart = 0;
		a[i].txEnd = 0;
		a[i].cdsStart = 0;
		a[i].cdsEnd = 0;
		a[i].exonCount = 0;
		a[i].ortho_id = -1;
		a[i].type = UNDEF;
		strcpy(a[i].gname, "");
		strcpy(a[i].sname, "");
		strcpy(a[i].chrname, "");
		strcpy(a[i].info, "");
	}
}

void initialize_blast_list(struct blast_list *a, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) {
		a[i].reg = assign_I(0, 1); 
		a[i].pid1 = (float) -1;
		a[i].pid2 = (float) -1;
		strcpy(a[i].name1, "");
		strcpy(a[i].name2, "");
		strcpy(a[i].info1, "");
		strcpy(a[i].info2, "");
	}
}

void initialize_exons(struct exons_list *a, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) {
		a[i].reg = assign_I(0, 1); 
		a[i].fid = -1;
		strcpy(a[i].name, "");
		strcpy(a[i].chrname, "");
	}
}

void quick_sort_dec_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list *h;
	float x;
	
	h = (struct g_list *) ckalloc(sizeof(struct g_list));
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
			assign_genes(h, 0, a[i]);
			assign_genes(a, i, a[j]);
			assign_genes(a, j, h[0]);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_genes(a, lo, j, mode);
	if (i < hi) quick_sort_dec_genes(a, i, hi, mode);
	free(h);
}

void quick_sort_inc_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list *h;
	float x;
	
	h = (struct g_list *) ckalloc(sizeof(struct g_list));
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
			assign_genes(h, 0, a[i]);
			assign_genes(a, i, a[j]);
			assign_genes(a, j, h[0]);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_genes(a, lo, j, mode);
	if (i < hi) quick_sort_inc_genes(a, i, hi, mode);
	free(h);
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

void assign_genes(struct g_list *new_genes, int loc,  struct g_list a)
{
  new_genes[loc].gid = a.gid;
  new_genes[loc].sid = a.sid;
  new_genes[loc].strand = a.strand;
  new_genes[loc].txStart = a.txStart;
  new_genes[loc].txEnd = a.txEnd;
  new_genes[loc].exonCount = a.exonCount;
  new_genes[loc].cdsStart = a.cdsStart;
  new_genes[loc].cdsEnd = a.cdsEnd;
	strcpy(new_genes[loc].gname, a.gname);
	strcpy(new_genes[loc].sname, a.sname);
	strcpy(new_genes[loc].chrname, a.chrname);
	strcpy(new_genes[loc].info, a.info);
  new_genes[loc].type = a.type;
  new_genes[loc].ortho_id = a.ortho_id;
}

void assign_exons(struct exons_list *new_exons, int loc, struct exons_list a)
{
  new_exons[loc].reg = assign_I(a.reg.lower, a.reg.upper);
  new_exons[loc].fid = a.fid;
	strcpy(new_exons[loc].name, a.name);
	strcpy(new_exons[loc].chrname, a.chrname);
}

void assign_blast_list(struct blast_list *new_a, int loc, struct blast_list a)
{
  new_a[loc].reg = assign_I(a.reg.lower, a.reg.upper);
  new_a[loc].pid1 = a.pid1;
  new_a[loc].pid2 = a.pid2;
	strcpy(new_a[loc].name1, a.name1);
	strcpy(new_a[loc].name2, a.name2);
	strcpy(new_a[loc].info1, a.info1);
	strcpy(new_a[loc].info2, a.info2);
}

void replace_gene_in_list(struct g_list *genes, int cur_id, int tmp_id)
{
  genes[cur_id].gid = genes[tmp_id].gid;
  genes[cur_id].sid = genes[tmp_id].sid;
  genes[cur_id].strand = genes[tmp_id].strand;
  genes[cur_id].txStart = genes[tmp_id].txStart;
  genes[cur_id].txEnd = genes[tmp_id].txEnd;
  genes[cur_id].exonCount = genes[tmp_id].exonCount;
  genes[cur_id].cdsStart = genes[tmp_id].cdsStart;
  genes[cur_id].cdsEnd = genes[tmp_id].cdsEnd;
	strcpy(genes[cur_id].gname, genes[tmp_id].gname);
	strcpy(genes[cur_id].sname, genes[tmp_id].sname);
	strcpy(genes[cur_id].chrname, genes[tmp_id].chrname);
	strcpy(genes[cur_id].info, genes[tmp_id].info);
  genes[cur_id].type = genes[tmp_id].type;
  genes[cur_id].ortho_id = genes[tmp_id].ortho_id;
}

int input_genes(FILE *f, struct g_list *genes, struct exons_list *exons)
{
	int i = -1, j = 0;
	char buf[10000], name[100], scf_name[100], chr_name[100];
	int num_genes = 0;
	int cur_exons_count = 0;
	int b = 0, e = 0;
	int old_b = 0, old_e = 0;

	strcpy(chr_name, "");
	strcpy(name, "");
	strcpy(buf, "");
	strcpy(scf_name, "");

	while(fgets(buf, 10000, f))
	{
		if( buf[0] == '#' ) {}
		else if( (buf[0] == '>') || (buf[0] == '<') )
		{
			if( i >= 0 ) {
				if( cur_exons_count == 0 ) {
					genes[i].cdsStart = j;
					exons[j].fid = i;
					exons[j].reg = assign_I(old_b, old_e);
					cur_exons_count++;
					j++;
				}
				genes[i].exonCount = cur_exons_count;
				genes[i].cdsEnd = j-1;
			}
			i++;
			cur_exons_count = 0;
			strcpy(chr_name, "");
			strcpy(name, "");
			if(buf[0] == '>') genes[i].strand = '+';
			else if(buf[0] == '<' ) genes[i].strand = '-';
			else fatalf("unexpected strand %c\n", buf[0]);

			if( strstr(buf, "(complement)") != NULL ) {
				if( sscanf(buf, "%*s %d %d %s %s %*s %s %*s", &b, &e, name, scf_name, chr_name) == 5 ) {
					strcpy(genes[i].chrname, chr_name);
					strcpy(genes[i].sname, scf_name);
				}
				else if( sscanf(buf, "%*s %d %d %s %s %*s", &b, &e, name, scf_name) == 4 ) {
					strcpy(genes[i].sname, scf_name);
				}
				else if( sscanf(buf,  "%*s %d %d %s %*s", &b, &e, name) != 3 ) {
					printf("wrong format in %s\n", buf);	
				}
				else strcpy(genes[i].sname, "");
			}
			else {
				if( sscanf(buf, "%*s %d %d %s %s %s %*s", &b, &e, name, scf_name, chr_name) == 5 ) {
					strcpy(genes[i].chrname, chr_name);
					strcpy(genes[i].sname, scf_name);
				}
				else if( sscanf(buf, "%*s %d %d %s %s %*s", &b, &e, name, scf_name) == 4 ) {
					strcpy(genes[i].sname, scf_name);
				}
				else if( sscanf(buf,  "%*s %d %d %s %*s", &b, &e, name) != 3 ) {
					printf("wrong format in %s\n", buf);	
				}
				else strcpy(genes[i].sname, "");
			}

			genes[i].type = ORF;
			genes[i].gid = i;
			genes[i].txStart = b;
			genes[i].txEnd = e;
			strcpy(genes[i].gname, name);
			old_b = b;
			old_e = e;
		}
		else {
			sscanf(buf, "%d %d", &b, &e);
			if( cur_exons_count == 0 ) genes[i].cdsStart = j;
			exons[j].fid = i;
			exons[j].reg = assign_I(b, e);
			cur_exons_count++;
			j++;
		}
	}

	if( i >= 0 ) { 
		if( cur_exons_count == 0 ) {
			genes[i].cdsStart = j;
			exons[j].fid = i;
			exons[j].reg = assign_I(old_b, old_e);
			cur_exons_count++;
			j++;
		}
		genes[i].exonCount = cur_exons_count;
		genes[i].cdsEnd = j-1;
	}

	i++;
	num_genes = i;

	return(num_genes);
}

int input_genes_blastx(FILE *f, struct blast_list *genes, int type, int cutoff_value, struct g_list *gff_genes, int num_gff_genes, struct exons_list *exons)
{
	char buf[10000], field[LEN_NAME];
	bool is_done = false;
	int b = 0, e = 0;
	int cur_b = 0, cur_e = 0;
	int len = 0, cur_len = 0, cds_len = 0;
	int num_hits = 0;
	int count = 0;
	float pid = (float) -1, sec_pid = (float) -1, cur_pid = (float) -1;
	char name[LEN_NAME], sec_name[LEN_NAME];
	char chr_name[LEN_NAME], sec_chr_name[LEN_NAME];
	char cur_name[LEN_NAME], cur_chr_name[LEN_NAME];
	int i = 0, j = 0;
	int num_genes = 0;
	int num_match = 0;
	struct I cur_reg, cds_reg, reg;
	int sid = 0, eid = 0;

	cur_reg = assign_I(0, 1);
	cds_reg = assign_I(0, 1);
	reg = assign_I(0, 1);
//	struct exons_list matches[2]; // name: gene_name, chr_name: chromosome name, reg: start and end, pid: percentage identity, value: e-value

	strcpy(field, "");
	strcpy(buf, "");

	while( (is_done == false) && fgets(buf, 10000, f) ) {
		while( (is_done == false) && (buf[0] == '#') ) {
			if( strstr(buf, "Query") != NULL ) {
				if( sscanf(buf, "%*s %*s %*s %s", field) != 1 ) {
					fatalf("wrong format in %s", buf);
				}
				else {
					split_b_and_e(field, &b, &e);	
					reg = assign_I(b, e);
					j = 0;
			    cur_reg = assign_I(gff_genes[j].txStart, gff_genes[j].txEnd);
			    while((j < num_gff_genes) && (equal(cur_reg, reg) == false)) 	
					{
      			j++;
						if( j < num_gff_genes ) cur_reg = assign_I(gff_genes[j].txStart, gff_genes[j].txEnd);
    			}

					if( j >= num_gff_genes ) {
						cds_reg =  assign_I(b, e);
					}
					else {
			      sid = gff_genes[j].cdsStart;
      			eid = gff_genes[j].cdsEnd;
      			if( exons[sid].reg.lower < exons[eid].reg.upper ) {
			        cds_reg = assign_I(exons[sid].reg.lower, exons[eid].reg.upper);
			      }
			      else {
			        if( gff_genes[j].strand == '-' )
 			      	{
			          cds_reg = assign_I(exons[eid].reg.lower, exons[sid].reg.upper);
			        }
			        else {
 			         fatalf("check CDS order in list: %d-%d\n", gff_genes[j].txStart, gff_genes[j].txEnd);
							}
		        }
		      }

					len = e - b;
					cds_len = width(cds_reg);
					num_match = 0;
					pid = (float) -1;
					sec_pid = (float) -1;
					strcpy(name, "");	
					strcpy(sec_name, "");	
					strcpy(chr_name, "");	
					strcpy(sec_chr_name, "");	
				}
			}
			else if( strstr(buf, "hit") != NULL ) {
				if( sscanf(buf, "%*s %d %*s", &num_hits) != 1 ) {
					fatalf("wrong format in %s", buf);
				}
			}
			if( !fgets(buf, 10000, f) ) is_done = true;
		}

		count = 0;
		while((is_done == false) && (buf[0] != '#')) {
			if( sscanf(buf, "%s %f %*s", cur_name, &cur_pid) != 2 ) {
				fatalf("wrong format in %s", buf);
			}

			if( cur_pid > (float)cutoff_value ) {
				parse_blastx_line(buf, cur_chr_name, &cur_b, &cur_e, type);
				cur_len = cur_e - cur_b;
				if( (((int)(((float)cur_len)/((float)len)*((float)100)) > BLAST_LEN_TH) && ((int)(((float)len)/((float)cur_len)*((float)100)) > BLAST_LEN_TH)) || (((int)(((float)cur_len)/((float)cds_len)*((float)100)) > BLAST_LEN_TH) && ((int)(((float)cds_len)/((float)cur_len)*((float)100)) > BLAST_LEN_TH) ) ) 
				{
					if( num_match  == 0 ) {
						pid = cur_pid;
						strcpy(name, cur_name);
						strcpy(chr_name, cur_chr_name);
						num_match++;
					}
					else if( num_match == 1 ) 
					{
						if( strcmp(name, cur_name) == 0 ) {}
						else if( cur_pid > pid ) 
						{
							sec_pid = pid;
							strcpy(sec_name, name);
							strcpy(sec_chr_name, chr_name);
							pid = cur_pid;
							strcpy(name, cur_name);
							strcpy(chr_name, cur_chr_name);
							num_match++;
						}
						else {
							sec_pid = cur_pid;
							strcpy(sec_name, cur_name);
							strcpy(sec_chr_name, cur_chr_name);
							num_match++;
						}
					}
					else {
						if( strcmp(name, cur_name) == 0 ) {}
						else if( strcmp(sec_name, cur_name) == 0 ) {}
						else if( cur_pid > pid ) {
							sec_pid = pid;
							strcpy(sec_name, name);
							strcpy(sec_chr_name, chr_name);
							pid = cur_pid;
							strcpy(name, cur_name);
							strcpy(chr_name, cur_chr_name);
						}
						else if( cur_pid > sec_pid ) {
							sec_pid = cur_pid;
							strcpy(sec_name, cur_name);
							strcpy(sec_chr_name, cur_chr_name);
						}
					}
				}
			}
			if( !fgets(buf, 10000, f) )is_done = true;
			count++;
		}
		
		if( (count != 0) && (num_hits > 0) && (num_match > 0) ) {
			genes[i].reg = assign_I(b, e);
			genes[i].pid1 = pid;
			strcpy(genes[i].name1, name);
			strcpy(genes[i].info1, chr_name);
			if( num_match >=2 ) {
				genes[i].pid2 = sec_pid;
				strcpy(genes[i].name2, sec_name);
				strcpy(genes[i].info2, sec_chr_name);
			}
			i++;
		}
	}

	num_genes = i;
	return(num_genes);
}

int input_genes_in_gff(FILE *f, struct g_list *genes, struct exons_list *exons)
{
	int i = -1, j = 0;
	char buf[10000], type[100], scf_name[100], chr_name[100];
	int num_genes = 0;
	int b = 0, e = 0;
	int old_b = 0, old_e = 0;
	char strand = '+';
	char gname[10000];
	int cur_exons_count = 0;
	char source[100];
	char column[10000];
	char item1[500], item2[500];
	char cur_gname[500], cds_gname[500];
	char name[500];

	strcpy(buf, "");
	strcpy(type, "");
	strcpy(scf_name, "");
	strcpy(chr_name, "");
	strcpy(source, "");
	strcpy(column, "");
	strcpy(item1, "");
	strcpy(item2, "");
	strcpy(cur_gname, "");
	strcpy(cds_gname, "");

  while( fgets(buf, 10000, f)) {
		strcpy(type, "");
		strcpy(scf_name, "");
		strcpy(chr_name, "");
		strcpy(source, "");
		strcpy(gname, "");
		strcpy(name, "");
		strcpy(column, "");
    if((buf[0] == '#') || (buf[0] == '-')) {
    }
    else {
     	if( (sscanf(buf, "%s %s %s %d %d %*s %c %*s %s %s %*s", scf_name, source, type, &b, &e, &strand, gname, chr_name) != 8) && ( sscanf(buf, "%s %s %s %d %d %*s %c %*s %s %*s", scf_name, source, type, &b, &e, &strand, gname) != 7 )) {
        fatalf("wrong gff format: %s", buf);
      }
      else if( (strcmp(type, "gene") == 0) || (strcmp(type, "match") == 0) || (strcmp(type, "match_part") == 0) ) {

				if( (num_genes > 0) && (cur_exons_count == 0) ) {
					genes[i].cdsStart = j;
					exons[j].fid = i;
					exons[j].reg = assign_I(old_b, old_e);
					cur_exons_count++;
					j++;
				} 

				if( i >= 0 ) {
					genes[i].exonCount = cur_exons_count;
					genes[i].cdsEnd = j-1;
				}

				i++;
				num_genes++;
				if( strstr(chr_name, "chr") != NULL ) {
					strcpy(genes[i].chrname, chr_name);
				}
				else {
					strcpy(genes[i].chrname, "");
				}

				strcpy(genes[i].sname, scf_name);
				if( strcmp(source, "SGD") == 0 ) {
					strcpy(column, gname);
          if( (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%s", item1, item2) != 2) ) {
            fatalf("wrong SGD gff column: %s", column);
          }
          else
          {
            if( sscanf(item1, "%*[^=]=%s", name) != 1 ) {
              fatalf("wrong gff column: %s", column);
            }
            strcpy(gname, name);
						strcpy(genes[i].info, item2);
						strcpy(cur_gname, gname);
          }
					strcpy(chr_name, "");
				}
				else {
					strcpy(cur_gname, "UNDEF");
				}

				strcpy(genes[i].gname, gname);
				genes[i].strand = strand;
				genes[i].gid = i;
				genes[i].txStart = b;
				genes[i].txEnd = e;
      	if( strcmp(type, "gene") == 0 ) {
					genes[i].type = ORF;
				}
				else if( strcmp(type, "match") == 0 ) {
					genes[i].type = MATCH;
				}
				else {
					genes[i].type = PARTIAL;
				}
				cur_exons_count = 0;
				old_b = b;
				old_e = e;
      }
      else if( strcmp(type, "CDS") == 0 ) {
				if( strcmp(source, "SGD") == 0 ) {
					strcpy(column, gname);

          if( (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%s", item1, item2) != 2) ) {
        		fatalf("wrong gff column: %s", column);
					}
					else {
						if( sscanf(item1, "%*[^=]=%s", name) != 1 ) {
        			fatalf("wrong gff column: %s", column);
						}
						else {
							if( sscanf(name, "%[^_]_%*s", item1) != 1 ) {
        				fatalf("wrong gff column: %s", name);
							}
						}
            strcpy(cds_gname, item1);
						strcpy(genes[i].info, item2);
					}
					strcpy(chr_name, "");
				}
				else strcpy(cds_gname, "UNDEF");

				if( strcmp(cur_gname, cds_gname) == 0 ) {
					if( cur_exons_count == 0 ) genes[i].cdsStart = j;
					exons[j].fid = i;
					exons[j].reg = assign_I(b, e);
					cur_exons_count++;
					j++;
				}
      }
    }
  }

	if( (num_genes > 0) && (cur_exons_count == 0) ) {
//		if( strcmp(cur_gname, cds_gname) == 0 ) {
			genes[i].cdsStart = j;
			exons[j].fid = i;
			exons[j].reg = assign_I(old_b, old_e);
			cur_exons_count++;
			j++;
//		}
	} 
	
	if( i >= 0 ) {
		genes[i].exonCount = cur_exons_count;
		genes[i].cdsEnd = j-1;
	}
	
	num_genes = i+1;
	return(num_genes);
}

int count_genes_blastx(FILE *f)
{
	char buf[10000];
	int num_genes = 0;
	int num_hits = 0;

	while( fgets(buf, 10000, f) ) {
		if( (buf[0] == '#') && (strstr(buf, "hits found") != NULL ) ) {
			if( sscanf(buf, "%*s %d %*s", &num_hits) != 1 ) {
				fatalf("wrong blastx output in %s", buf);
			}
			else if( num_hits > 0 ) {
				num_genes++;
			}
		}
	}
	return(num_genes);
}

int count_genes(FILE *f, int *num_exons)
{
	char buf[10000];
	int cur_exons_count = 0;
	int num_genes = 0;

	while( fgets(buf, 10000, f) ) {
		if( buf[0] == '#' ) {}
		else if( ( buf[0] == '>' ) || ( buf[0] == '<' ) ) {
			if( cur_exons_count == 0 ) {
          *num_exons = *num_exons + 1;
			}
			num_genes++;
			cur_exons_count = 0;
		}
		else {
			*num_exons = *num_exons + 1;
			cur_exons_count++;
		}
	}

	if( num_genes >= 0 ) { 
		if( cur_exons_count == 0 ) {
			*num_exons = *num_exons + 1;
		}
	}

	return(num_genes);
}

int count_genes_in_gff(FILE *f, int *num_exons)
{
	char buf[10000];
	int cur_exons_count = 0;
	int num_genes = 0;
	int b = 0, e = 0;
	char type[100];

  while( fgets(buf, 10000, f)) {
    if((buf[0] == '#') || (buf[0] == '-')) {
    }
    else {
      if( sscanf(buf, "%*s %*s %s %d %d %*s", type, &b, &e) != 3 ) {
        fatalf("wrong gff format: %s", buf);
      }
      else if( (strcmp(type, "gene") == 0) || (strcmp(type, "match") == 0) || (strcmp(type, "match_part") == 0) ) {
				if( (num_genes > 0) && (cur_exons_count == 0) ) {
					*num_exons = *num_exons + 1;
				}
				num_genes++;
				cur_exons_count = 0;
			}
      else if( strcmp(type, "CDS") == 0 ) {
				cur_exons_count++;
				*num_exons = *num_exons + 1;
			}
		}
	}

  if( (num_genes > 0) && (cur_exons_count == 0) ) {
    *num_exons = *num_exons + 1;
  }

	return(num_genes);
}

int rm_redun_genes(struct g_list *genes, int from, int to)
{
	int i = 0, j = 0;
	struct I reg, tmp_reg;

	reg = assign_I(0, 1);
	tmp_reg = assign_I(0, 1);

	for( i = from; i < to; i++ ) {
		if( genes[i].type == REDUN ) {}
		else {
			reg = assign_I(genes[i].txStart, genes[i].txEnd);
			j = i+1;	
			tmp_reg = assign_I(genes[j].txStart, genes[j].txEnd);
			while( (j <= to) && (overlap(reg, tmp_reg) == true) ) {
				if( genes[j].type == REDUN ) {}
				else if(near_equal(reg, tmp_reg, OV_TH) == true) {
					if( (genes[i].type == ORF) && (genes[j].type == ORF ) ) {
						if( equal(reg, tmp_reg) == true ) {
							if( genes[i].exonCount < genes[j].exonCount ) {
								genes[j].type = REDUN;
							}
							else if (genes[i].exonCount > genes[j].exonCount) {
								genes[i].type = REDUN;
							}
							else if( strstr(genes[i].gname, "UNDEF") != NULL ) {
								genes[i].type = REDUN;
							}
							else {
								genes[j].type = REDUN;
							}
						}
					}
					else if( genes[i].type == ORF ) {
						genes[j].type = REDUN;
					}
					else if( genes[j].type == ORF ) {
						genes[i].type = REDUN;
					}
					else {
						if( genes[i].exonCount < genes[j].exonCount ) {
							genes[j].type = REDUN;
						}
						else if (genes[i].exonCount > genes[j].exonCount) {
							genes[i].type = REDUN;
						}
						else genes[i].type = REDUN;
					}
				}
				else if( subset(tmp_reg, reg) == true ) {
					if( genes[i].type == PARTIAL ) {
						genes[i].type = REDUN;
					}
				}
	
				j++;	
				if( j <= to ) tmp_reg = assign_I(genes[j].txStart, genes[j].txEnd);
			}
		}
	}

	j = 0;
	for( i = from; i <= to; i++ ) {
		if( genes[i].type != REDUN ) {
			if( j != i ) replace_gene_in_list(genes, j, i);
			j++;
		}
	}

	return(j);
}

int count_exons(struct g_list *genes, int from, int to)
{
	int i = 0;
	int count = 0;

	for( i = from; i <=to; i++ ) {
		count = count + genes[i].cdsEnd - genes[i].cdsStart + 1;
	}
	
	return(count);
}

void write_in_gff(struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons)
{
	int i = 0, j = 0;
	char chrname[LEN_NAME], gname[LEN_NAME];

	strcpy(chrname, "");
	strcpy(gname, "");

	for( i = 0; i < num_genes; i++ ) {
		if( (strcmp(genes[i].chrname, "") != 0) && (strcmp(genes[i].gname, "") != 0) ) {
			strcpy(chrname, genes[i].chrname);
			strcpy(gname, genes[i].gname);
			if( strstr(chrname, gname) != NULL ) {}
			else {
				strcpy(genes[i].chrname, gname);
				strcat(genes[i].chrname, ",SacCer,");
				strcat(genes[i].chrname, chrname);
				strcat(genes[i].chrname, ",LOW_MATCH");
			}
		}

		if( genes[i].type == ORF ) {
			if( genes[i].ortho_id == DOUBLE_ORFS ) {
				if( strcmp(genes[i].info, "") == 0 ) {
					fatalf("%s: info field supposed not empty\n", genes[i].gname);
				}
				else {
					printf("%s maker gene %d %d . %c . %s;%s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname, genes[i].info);
				}
			}
			else if( strcmp(genes[i].chrname, "") != 0) {
				printf("%s maker gene %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname);
			}
			else if( strcmp(genes[i].gname, "") != 0) {
				printf("%s maker gene %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].gname);
			}
			else {
				printf("%s maker gene %d %d . %c . UNDEF\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand);
			}

			for( j = genes[i].cdsStart; j <= genes[i].cdsEnd; j++ ) {
				if( j >= num_exons ) {
					fatalf("%d excceeds range [0,%d]\n", j, num_exons);
				}

				if( genes[i].ortho_id == DOUBLE_ORFS ) {
					if( strcmp(genes[i].info, "") == 0 ) {
						fatalf("%s: info field supposed not empty\n", genes[i].gname);
					}
					else {
						printf("%s maker CDS %d %d . %c . %s;%s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper, genes[i].strand, genes[i].chrname, genes[i].info);
					}
				}
				else if( strcmp(genes[i].chrname, "") != 0 ) {
					printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper, genes[i].strand, genes[i].chrname);
				}
				else if( strcmp(genes[i].gname, "") != 0) {
					printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper, genes[i].strand, genes[i].gname);
				}
				else {
					printf("%s maker CDS %d %d . %c . UNDEF\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper, genes[i].strand);
				}
			}
		}
		else if( genes[i].type == MATCH ) {
			if( genes[i].ortho_id == DOUBLE_ORFS ) {
				if( strcmp(genes[i].info, "") == 0 ) {
					fatalf("%s: info field supposed not empty\n", genes[i].gname);
				}
				else {
					printf("%s maker match %d %d . %c . %s;%s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname, genes[i].info);
				}
			}
			else if( strcmp(genes[i].chrname, "") != 0 ) {
				printf("%s maker match %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname);
			}
			else {
				printf("%s maker match %d %d . %c . UNDEF\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand);
			}
		}
		else if( genes[i].type == PARTIAL ) {
			if( genes[i].ortho_id == DOUBLE_ORFS ) {
				if( strcmp(genes[i].info, "") == 0 ) {
					fatalf("%s: info field supposed not empty\n", genes[i].gname);
				}
				else {
					printf("%s maker match_part %d %d . %c . %s;%s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname, genes[i].info);
				}
			}
			else if( strcmp(genes[i].chrname, "") != 0 ) {
				printf("%s maker match_part %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, genes[i].chrname);
			}
			else {
				printf("%s maker match_part %d %d . %c . UNDEF\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand);
			}
		}
		else {
//			fatalf("%s:%d-%d type:%d undefined\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].type);
		}
	}
}

int rm_overlap_genes(struct g_list *genes, struct exons_list *exons, int from, int to)
{
	int i = 0, j = 0;
	struct I reg, tmp_reg;
	struct I exon_reg, exon_tmp;
	double evalue1 = (double) 1, evalue2 = (double) 1;
	float pid1 = (float) 0, pid2 = (float) 0;
	int num_values = 0;
	char res1[LEN_NAME], res2[LEN_NAME];
	bool is_first = true;
	int len1 = 0, len2 = 0;
	int sid = 0, eid = 0;

	strcpy(res1, "");
	strcpy(res2, "");

	reg = assign_I(0, 1);
	tmp_reg = assign_I(0, 1);
	exon_reg = assign_I(0, 1);
	exon_tmp = assign_I(0, 1);

	for( i = from; i < to; i++ ) {
		if( strstr(genes[i].gname, "LOW_MATCH") != NULL ) {
			strcpy(genes[i].gname, "UNDEF");
		}
	}

	for( i = from; i < to; i++ ) {
		if( genes[i].type == REDUN ) {}
		else {
			reg = assign_I(genes[i].txStart, genes[i].txEnd);
			sid = genes[i].cdsStart;
			eid = genes[i].cdsEnd;
			if( exons[sid].reg.lower < exons[eid].reg.upper ) {
				exon_reg = assign_I(exons[sid].reg.lower, exons[eid].reg.upper);	
			}
			else {
				if( genes[i].strand == '-' )
				{
					exon_reg = assign_I(exons[eid].reg.lower, exons[sid].reg.upper);
				}
				else {
					fatalf("check CDS order in list: %d-%d\n", genes[i].txStart, genes[i].txEnd);
				}
			}

			is_first = true;
			if( strcmp(genes[i].gname, "UNDEF") != 0 ) {
				num_values = get_column_in_gff_nth_column(genes[i].gname, 8, res1, res2); // 8th column
				pid1 = get_pid(res1, res2, num_values); 
				if( pid1 >= PID_CUTOFF ) {
					num_values = get_column_in_gff_nth_column(genes[i].gname, 7, res1, res2); // 7th column
					if( (num_values == 0) || (num_values > 2) ) {
						fatalf("expected not UNDEF in %s\n", genes[i].gname);
					}
					else {
						evalue1 = get_evalue(res1, res2, num_values, &is_first); 
						if( evalue1 <= MIN_EVAL_CUTOFF ) evalue1 = (double) 0;
						num_values = get_column_in_gff_nth_column(genes[i].gname, 8, res1, res2); // 8th column
						if( is_first == true ) {
							pid1 = atof(res1);
						}
						else pid1 = atof(res2);
						len1 = get_len_diff(genes[i].gname, exon_reg, is_first);
					}
	
					if( evalue1 > EVAL_CUTOFF ) {
						if( (width(exon_reg) < MIN_NU_LEN) && ((float)len1/(float)(width(exon_reg)) < MIN_LEN_RATIO) && ((int)pid1 > HIGH_PID) ) {
						}
						else {
							genes[i].type = REDUN;
						}
					}
				}
				else {
					genes[i].type = REDUN;
				}
			}
			j = i+1;	
			tmp_reg = assign_I(genes[j].txStart, genes[j].txEnd);
			sid = genes[j].cdsStart;
			eid = genes[j].cdsEnd;
			if( exons[sid].reg.lower < exons[eid].reg.upper ) {
				exon_tmp = assign_I(exons[sid].reg.lower, exons[eid].reg.upper);	
			}
			else {
				if( genes[j].strand == '-' )
				{
					exon_tmp = assign_I(exons[eid].reg.lower, exons[sid].reg.upper);
				}
				else {
					fatalf("check CDS order in list: %d-%d\n", genes[j].txStart, genes[j].txEnd);
				}
			}
			while( (j <= to) && (proper_overlap(reg, tmp_reg) == true) ) {
				if( genes[j].type == REDUN ) {}
				else if(mostly_overlap(reg, tmp_reg, CUT_RATIO) == true) {
					if( (genes[i].type == ORF ) && (genes[j].type != ORF) ) {
						genes[j].type = REDUN;	
					}
					else if( (genes[i].type != ORF ) && (genes[j].type == ORF) ) {
						genes[i].type = REDUN;	
					}
					else if( (strcmp(genes[i].gname, "UNDEF") == 0 ) && (strcmp(genes[j].gname, "UNDEF") == 0 ) ) {
						if( width(tmp_reg) >= width(reg) ) genes[i].type = REDUN;
						else genes[j].type = REDUN;
					}
					else if( strcmp(genes[i].gname, "UNDEF") == 0 ) {
						genes[i].type = REDUN;
					}
					else if( strcmp(genes[j].gname, "UNDEF") == 0 ) {
						genes[j].type = REDUN;
					}
					else {
						is_first = true;
						num_values = get_column_in_gff_nth_column(genes[j].gname, 8, res1, res2); // 8th column
						pid2 = get_pid(res1, res2, num_values); 
						if( pid2 >= PID_CUTOFF ) {
							num_values = get_column_in_gff_nth_column(genes[j].gname, 7, res1, res2);
							evalue2 = get_evalue(res1, res2, num_values, &is_first); 
							if( evalue2 <= MIN_EVAL_CUTOFF ) evalue2 = (double) 0;
							num_values = get_column_in_gff_nth_column(genes[j].gname, 8, res1, res2); // 8th column
							if( is_first == true ) {
								pid2 = atof(res1);
							}
							else pid2 = atof(res2);
							len2 = get_len_diff(genes[j].gname, exon_tmp, is_first);

							if( evalue2 > EVAL_CUTOFF ) {
								if( (width(exon_tmp) < MIN_NU_LEN) && ((float)len2/(float)(width(exon_tmp)) < MIN_LEN_RATIO) && ((int)pid2 > HIGH_PID) ) {
								}
								else {
									genes[j].type = REDUN;
								}
							}
							else {
								if( evalue1 > evalue2 ) genes[i].type = REDUN;
								else if( evalue1 < evalue2 ) genes[j].type = REDUN;
								else {
									if( len1 < len2 ) genes[j].type = REDUN;
									else if( len1 > len2 ) genes[i].type = REDUN;
									else {
										if( pid1 < pid2 ) genes[i].type = REDUN;
										else if( pid1 > pid2) genes[j].type = REDUN;
										else {
											if( width(tmp_reg) >= width(reg) ) genes[i].type = REDUN;
											else genes[j].type = REDUN;
										}
									}
								}
							}
						}
						else {
							genes[j].type = REDUN;
						}
					}
				}

				j++;	
				if( j <= to ) {
					tmp_reg = assign_I(genes[j].txStart, genes[j].txEnd);
					sid = genes[j].cdsStart;
					eid = genes[j].cdsEnd;
					if( exons[sid].reg.lower < exons[eid].reg.upper ) {
						exon_tmp = assign_I(exons[sid].reg.lower, exons[eid].reg.upper);	
					}
					else {
						if( genes[j].strand == '-' )
						{
							exon_tmp = assign_I(exons[eid].reg.lower, exons[sid].reg.upper);
						}
						else {
							fatalf("check CDS order in list: %d-%d\n", genes[j].txStart, genes[j].txEnd);
						}
					}
				}
			}
		}
	}

	j = 0;
	for( i = from; i <= to; i++ ) {
		if( genes[i].type != REDUN ) {
			if( j != i ) replace_gene_in_list(genes, j, i);
			j++;
		}
	}

	return(j);
}

int get_column_in_gff_nth_column(char *column, int nth, char *res1, char *res2)
{
	int i = 0, j = 0;
	char value[LEN_NAME];
	int num_columns = 0;
	int num_vals = 0;

	strcpy(value, "");	
	strcpy(res1, "");
	strcpy(res2, "");
	while( (column[i] != '\0') && (column[i] != '\n') ) {
		if( column[i] == ',' ) {
			value[j] = '\0';
			j = 0;
			num_columns++;
			if( num_columns == nth ) {
				if( num_vals == 0 ) {
					strcpy(res1, value);
					num_vals++;
				}
				else if( num_vals == 1 ) {
					strcpy(res2, value);
					num_vals++;
				}
			}
		}
		else if( column[i] == ';' ) {
			value[j] = '\0';
			j = 0;
			num_columns++;
			if( num_columns == nth ) {
				if( num_vals == 0 ) {
					strcpy(res1, value);
					num_vals++;
				}
				else if( num_vals == 1 ) {
					strcpy(res2, value);
					num_vals++;
				}
			}
		}
		else {
			value[j] = column[i];
			j++;
		}

	
		if( column[i] == ';' ) {
			num_columns = 0;
		}
		i++;
	}
	num_columns++;
	value[j] = '\0';
	if( num_columns == nth ) {
		if( num_vals == 0 ) {
			strcpy(res1, value);
			num_vals++;
		}
		else if( num_vals == 1 ) {
			strcpy(res2, value);
			num_vals++;
		}
	}

	return(num_vals);
}

void write_orfs_in_other_splicing(struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons, int cut_off)
{
	int i = 0, j = 0;
	char gname[2][LEN_NAME];
	char seq_name[2][LEN_NAME];
	char chr[2][LEN_NAME];
	int b[2], e[2];
	float pid[2];
	char line[2][LEN_LONG_NAME];
	int max_pid = 0;
	int cur_b = 0, cur_e = 0;
	char cur_gname[LEN_NAME];
	int count = 0;
	 
	for( i = 0; i < 2; i++ ) {
		strcpy(gname[i], "");
		strcpy(seq_name[i], "");
		strcpy(chr[i], "");
		b[i] = 0;
		e[i] = 0;
		pid[i] = (float) 0;
	}

	for( i = 0; i < num_genes; i++ ) 
	{
		if( strstr(genes[i].gname, ";") != NULL ) {
			if( sscanf(genes[i].gname, "%[^;];%s", line[0], line[1]) != 2 ) 
			{
				fatalf("wrong format in the last column: %s\n", genes[i].gname);
			}	
		}
		else {
			strcpy(line[0], genes[i].gname);
		}
		
		if( strstr(genes[i].gname, "UNDEF") != NULL ) {
			max_pid = 0;
			strcpy(cur_gname, "UNDEF");
		}
		else {
			parse_last_column_gff(line[0], gname[0], seq_name[0], chr[0], &b[0], &e[0], &pid[0]);
			if( strstr(genes[i].gname, ";") != NULL ) {
				parse_last_column_gff(line[1], gname[1], seq_name[1], chr[1], &b[1], &e[1], &pid[1]);
			}

			if( pid[0] > pid[1] ) {
				max_pid = (int) pid[0];
				strcpy(cur_gname, gname[0]);
			}
			else {
				max_pid = (int) pid[1];
				strcpy(cur_gname, gname[1]);
			}
		}

		if( ((genes[i].cdsEnd - genes[i].cdsStart) >= 1) && (max_pid < cut_off)) {
			if( genes[i].strand == '+' ) {
				for( j = genes[i].cdsStart; j <= genes[i].cdsEnd; j++ ) {
					if( j >= num_exons ) {
						fatalf("invalid exon index %d beyond the number of exons %d\n", j, num_exons);
					}

					if( j == genes[i].cdsStart ) {
						cur_b = genes[i].txStart;
						cur_e = exons[j+1].reg.lower;	
					}
					else if( j == genes[i].cdsEnd ) {
						cur_b = exons[j-1].reg.upper;
						cur_e = genes[i].txEnd;
					}
					else {
						cur_b = exons[j-1].reg.upper;
						cur_e = exons[j+1].reg.lower;	
					}
					count++;
					printf("> %d %d %s-%d\n", cur_b, cur_e, genes[i].sname, count);		
					printf("%d %d\n", cur_b, cur_e);		
				}
			}
			else {
				for( j = genes[i].cdsEnd; j >= genes[i].cdsStart; j-- ) {
					if( j >= num_exons ) {
						fatalf("invalid exon index %d beyond the number of exons %d\n", j, num_exons);
					}

					if( j == genes[i].cdsStart ) {
						cur_b = exons[j+1].reg.upper;	
						cur_e = genes[i].txEnd;
					}
					else if( j == genes[i].cdsEnd ) {
						cur_b = genes[i].txStart;
						cur_e = exons[j-1].reg.lower;
					}
					else {
						cur_b = exons[j+1].reg.upper;
						cur_e = exons[j-1].reg.lower;	
					}
					count++;
					printf("> %d %d %s-%d\n", cur_b, cur_e, genes[i].sname, count);		
					printf("%d %d\n", cur_b, cur_e);		
				}
			}
		}  
	}
}

float get_pid(char *res1, char *res2, int num_values) 
{
	float val1 = (float) 0, val2 = (float) 0;

	if( num_values == 1 ) {
		val1 = atof(res1);
	}
	else if( num_values == 2) {
		val1 = atof(res1);
		val2 = atof(res2);
	}

	if( val1 > val2 ) return(val1);
	else return(val2);
}

double get_evalue(char *res1, char *res2, int num_values, bool *is_first) 
{
	double val1 = (double) 100, val2 = (double) 100;

	*is_first = true;
	if( num_values == 1 ) {
		val1 = AtoF(res1);
	}
	else if( num_values == 2) {
		val1 = AtoF(res1);
		val2 = AtoF(res2);
	}

	if( val1 <= val2 ) return(val1);
	else {
		*is_first = false;
		return(val2);
	}
}

int get_len_diff(char *gname, struct I reg, bool is_first)
{
	int num_values = 0;
	int b = 0, e = 0;
	struct I cur_reg;
	char res1[LEN_NAME], res2[LEN_NAME];
	int len = 0;

	cur_reg = assign_I(0, 1);
	num_values = get_column_in_gff_nth_column(gname, 4, res1, res2); 
	if( is_first == true ) {
		b = atoi(res1);
	}
	else b = atoi(res2);

	num_values = get_column_in_gff_nth_column(gname, 5, res1, res2); 
	if( is_first == true ) {
		e = atoi(res1);
	}
	else e = atoi(res2);

	cur_reg = assign_I(b, e);
	len = abs(width(reg) - width(cur_reg));

	return(len);
}

void print_item_gff(struct g_list gene, struct exons_list *exons)
{
	int i = 0;

	if( gene.type == ORF ) {
		printf("%s\tannot\tgene\t%d\t%d\t.\t%c\t.\t%s\n", gene.sname, gene.txStart, gene.txEnd, gene.strand, gene.gname);
		for( i = gene.cdsStart; i <= gene.cdsEnd; i++ ) {
			printf("%s\tannot\tCDS\t%d\t%d\t.\t%c\t.\t%s\n", gene.sname, exons[i].reg.lower, exons[i].reg.upper, gene.strand, gene.gname);
		}
	}
	else if( gene.type == MATCH ) {
		printf("%s\tannot\tmatch\t%d\t%d\t.\t%c\t.\t%s\n", gene.sname, gene.txStart, gene.txEnd, gene.strand, gene.gname);
	}
	else if( gene.type == PARTIAL ) {
		printf("%s\tannot\tmatch_part\t%d\t%d\t.\t%c\t.\t%s\n", gene.sname, gene.txStart, gene.txEnd, gene.strand, gene.gname);
	}
}
