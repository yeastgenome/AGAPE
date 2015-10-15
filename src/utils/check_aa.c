#include "main.h"
#include "util.h"
#include "tokens.h"
#include "util_genes.h"
#include "util_i.h"

char compl[128];

void print_genes(int b, int e, struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons, char *type, int shift_b, int shift_e, bool is_internal_stops);

int main(int argc, char **argv) {
	FILE *fp;
	char buf[1000], info[1000];
	struct g_list *genes;
	struct exons_list *exons;
	int num_genes = 0, num_exons = 0;
	int i = 0, j = 0;
	bool is_correct_stop = false;
	bool is_wrong_splicing = true;
	bool is_correct_start = false;
	int count = 0;
	int num_lines = 0;
	char last_aa = '\0';
	int len = 0;
	int b = 0, e = 0;
	int shift_b = 0, shift_e = 0;
	bool is_internal_stops = false;

	if( argc == 4 ) {
		if(strcmp(argv[3], "INTERNAL_STOPS") == 0) {
			is_internal_stops = true;	
		}
		else 
			fatalf("args: gff dna");
	}
	else if (argc != 3)
		fatalf("args: gff dna");
	compl['a'] = compl['A'] = 'T';
	compl['c'] = compl['C'] = 'G';
	compl['g'] = compl['G'] = 'C';
	compl['t'] = compl['T'] = 'A';

	strcpy(buf, "");

  if((fp = ckopen(argv[1], "r")) == NULL )
  {
    fatalf("Cannot open file %s\n", argv[2]);
  }
  else {
    num_genes = count_genes_in_gff(fp, &num_exons);
    if( num_genes > 0 ) {
      genes = (struct g_list *) ckalloc(num_genes * sizeof(struct g_list));
      if( num_exons < num_genes ) num_exons = num_genes;
      exons = (struct exons_list *) ckalloc(num_exons * sizeof(struct exons_list));
      initialize_genes(genes, num_genes);
      initialize_exons(exons, num_exons);
    }
  }

  fseek(fp, 0, SEEK_SET);
  count = input_genes_in_gff(fp, genes, exons);
  if( count != num_genes ) {
    fatalf("gene counter error in %s\n", argv[1]);
  }

  if( num_genes > 0 ) {
    quick_sort_inc_genes(genes, 0, num_genes-1, POS_BASE);
  }
  i = 0;
  while( i < num_genes ) {
    j = 0;
    while( ((i+j) < num_genes) && (genes[i].txStart == genes[i+j].txStart )) j++;
    quick_sort_dec_genes(genes, i, i+j-1, LEN_BASE);
    i = i+j;
  }
	
	fclose(fp);
	fp = ckopen(argv[2], "r");

	count = 0;
	while (fgets(buf, 1000, fp)) {
		if( buf[0] == '>' ) {
			if( count != 0 ) {
				if( (is_correct_start == true) && (is_wrong_splicing == false) && (is_correct_stop == true ) ) {
					if( is_internal_stops == false) 
						print_genes(b, e, genes, num_genes, exons, num_exons, "gene", shift_b, shift_e, is_internal_stops);
				}				
				else {
					if( is_internal_stops == false) 
						print_genes(b, e, genes, num_genes, exons, num_exons, "match", shift_b, shift_e, is_internal_stops);
					else 
						print_genes(b, e, genes, num_genes, exons, num_exons, "gene", shift_b, shift_e, is_internal_stops);
				}
			}

			if(sscanf(buf, "%*s %s %*s", info) != 1 ) {
			}
			else {
				split_b_and_e(info, &b, &e);
			}
			count++;
			num_lines = 0;
			shift_b = 0;
			shift_e = 0;
			is_correct_start = false;
			is_correct_stop = false;
			is_wrong_splicing = false;
		}
		else if( isspace(buf[0]) ) {}
		else {
			if( num_lines == 0 ) {
				if( buf[0] == 'M' ) is_correct_start = true;
				else {
					i = 1;
					while( (is_correct_start == false) && (buf[i] != '\0') && (buf[i] != '\n') ) {
						if( buf[i] == 'M' ) {
							is_correct_start = true;
							shift_b = i*3;
						}
						i++;
					}
				}
			}
			
			len = strlen(buf);
			if( isspace(buf[len-1]) ) {
				last_aa = buf[len-2];
				buf[len-2] = '\0';
				i = len-3;
			}
			else {
				last_aa = buf[len-1];
				buf[len-1] = '\0';
				i = len-2;
			}

			shift_e = 0;
			if( (strstr(buf, "*") != 0) || (last_aa == '*') ) {
				if( (strstr(buf, "*") != 0) && (last_aa == '*') ) {
					j = i;
					while( j >=0 ) {
						if( buf[j] == '*' ) {
							shift_e = (i-j+1) * 3;
						}		
						j--;
					}	
				}

				if( is_correct_stop == true ) {
					is_wrong_splicing = true;
				}
				is_correct_stop = true;
			}

			num_lines++;	
		}
	}
	
	if( count != 0 ) {
		if( (is_correct_start == true) && (is_wrong_splicing == false) && (is_correct_stop == true ) ) {
			if( is_internal_stops == false) 
				print_genes(b, e, genes, num_genes, exons, num_exons, "gene", shift_b, shift_e, is_internal_stops);
		}				
		else 
			if( is_internal_stops == false) 
				print_genes(b, e, genes, num_genes, exons, num_exons, "match", shift_b, shift_e, is_internal_stops);
			else 
				print_genes(b, e, genes, num_genes, exons, num_exons, "gene", shift_b, shift_e, is_internal_stops);
	}

	if( num_genes > 0 ) {
		free(genes);
		free(exons);
	}
	fclose(fp);
	return EXIT_SUCCESS;
}

void print_genes(int b, int e, struct g_list *genes, int num_genes, struct exons_list *exons, int num_exons, char *type, int shift_b, int shift_e, bool is_internal_stops)
{
	int cur_b = 0, cur_e = 0;
	bool is_same = false;
	int i = 0, j = 0;
	char cur_type[100];
	char gname[1000];
	char column[10000];
	char item1[500];
	char item2[500];
		
	strcpy(gname, "");
	strcpy(column, "");
	strcpy(item1, "");
	strcpy(item2, "");

	strcpy(cur_type, "");
	strcpy(cur_type, type);

	while( (is_same == false) && (i < num_genes) ) {
		cur_b = genes[i].txStart;
		cur_e = genes[i].txEnd;	

		if( (b == cur_b) && (e == cur_e) ) {

			if( strcmp(type, "gene") == 0 ) {
				j = genes[i].cdsStart;
				if( shift_b <= width(exons[j].reg) ) {}
				else {
					strcpy(cur_type, "match");
				}

				j = genes[i].cdsEnd;
				if( shift_e <= width(exons[j].reg) ) {}
				else {
					strcpy(cur_type, "match");
				}
			}

			if( is_internal_stops == true ) {
				strcpy(column, genes[i].gname);
				strcpy(gname, genes[i].gname);			
			}
			else {
				strcpy(gname, "UNDEF");
			}

			if( strcmp(cur_type, "gene") == 0 ) {
     		printf("%s maker gene %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, gname);
			}	
			else if( strcmp(cur_type, "match") == 0 ) {
     		printf("%s maker match %d %d . %c . %s\n", genes[i].sname, genes[i].txStart, genes[i].txEnd, genes[i].strand, gname);
			}
			else {
				fatalf("unknown type: %s\n", type);
			}
    }

		if( (b == cur_b) && (e == cur_e) && strcmp(type, "gene") == 0 ) {
   	 for( j = genes[i].cdsStart; j <= genes[i].cdsEnd; j++ ) {
				if( j >= num_exons ) {
					fatalf("%d excceeds range [0,%d]\n", j, num_exons);
				}	
				
				if( (j == genes[i].cdsStart) && (j == genes[i].cdsEnd) ) {
					if( genes[i].strand == '+' ) {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower+shift_b, exons[j].reg.upper-shift_e, genes[i].strand, gname);
					}
					else {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower+shift_e, exons[j].reg.upper-shift_b, genes[i].strand, gname);
					}
				}
				else if( j == genes[i].cdsStart ) {
					if( genes[i].strand == '+' ) {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower+shift_b, exons[j].reg.upper, genes[i].strand, gname);
					}
					else {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper-shift_b, genes[i].strand, gname);
					}
				}
				else if( j == genes[i].cdsEnd ) {
					if( genes[i].strand == '+' ) {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper-shift_e, genes[i].strand, gname);
					}
					else {
						printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower+shift_e, exons[j].reg.upper, genes[i].strand, gname);
					}
				}
				else {
					printf("%s maker CDS %d %d . %c . %s\n", genes[i].sname, exons[j].reg.lower, exons[j].reg.upper, genes[i].strand, gname);
				}
			}
		}
		
		i++;
	}
}
