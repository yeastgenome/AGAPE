// check_splicing_signal.c -- branchpoint
// starting position: 1 and both ends are included
#include "main.h"
#include "util_genes.h"
#include "util.h"
#include "util_i.h"
#include "seq.h"

#define NUM_BP_SEQ 9

char compl[128];
char **branchpoints;

bool check_strstr(char *cur_seq, char **branchpoints, int num_bp_seq);
bool check_introns(struct g_list *genes, int gid, struct exons_list *exons, char *seq, int seq_len);

int main(int argc, char **argv) {
	SEQ *sf;
	uchar *s;
	FILE *f;
	char chr_name[100], info[1000], dir;
	int i = 0, j = 0, k = 0, B = 0, E = 0;
	int max_len = 0;
	char *cur_seq;
	int seq_len = 0;
	bool is_correct_splicing = false;
	int num_genes = 0;
	int num_exons = 0;
	int num = 0;
	struct g_list *genes;
	struct exons_list *exons;
	bool no_branchpoint = false;

	if( argc == 5 ) {
		
	}	
	else if( argc == 4 ) {
		if( strcmp( argv[3], "NO_BRANCHPOINT") == 0 ) {
			no_branchpoint = true;
		}
		else {
			fatalf("args: fasta gff (NO_BRANCHPOINT)");
		}
	}
	else if (argc != 3)
		fatalf("args: fasta gff (NO_BRANCHPOINT)");

	if((f = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		num_genes = count_genes_in_gff(f, &num_exons);
		if( num_genes > 0 ) {
			genes = (struct g_list *) ckalloc(num_genes * sizeof(struct g_list));	
			if( num_exons < num_genes ) num_exons = num_genes;
			exons = (struct exons_list *) ckalloc(num_exons * sizeof(struct exons_list));	
			initialize_genes(genes, num_genes);
			initialize_exons(exons, num_exons);
		}
	}
	fseek(f, 0, SEEK_SET);
	
	branchpoints = (char **) ckalloc(sizeof(char *) * NUM_BP_SEQ);
	for( i = 0; i < NUM_BP_SEQ; i++ ) {
		branchpoints[i] = (char *) ckalloc(sizeof(char) * 8);	
	}
	
	strcpy(branchpoints[0], "AACTAAC");
	strcpy(branchpoints[1], "AATTAAC");
	strcpy(branchpoints[2], "CACTAAC");
	strcpy(branchpoints[3], "GACTAAC");
	strcpy(branchpoints[4], "TACTAAC");
	strcpy(branchpoints[5], "TACTAAT");
	strcpy(branchpoints[6], "TATTAAC");
	strcpy(branchpoints[7], "TGCTAAC");
	strcpy(branchpoints[8], "GATTAAC");

	num = input_genes_in_gff(f, genes, exons);	
	if( num != num_genes ) {
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
	fclose(f);

	compl['a'] = compl['A'] = 'T';
	compl['c'] = compl['C'] = 'G';
	compl['g'] = compl['G'] = 'C';
	compl['t'] = compl['T'] = 'A';
	sf = seq_get(argv[1]);
	s = SEQ_CHARS(sf) - 1;
	seq_len = SEQ_LEN(sf);

	for( i = 0; i < num_genes; i++ ) {
		B = genes[i].txStart;
		E = genes[i].txEnd;
		if( E > seq_len ) {
			fatalf("gene boundary [%d,%d] over the sequence length %d\n", B, E, seq_len);
		}

		if( (E - B + 1) > max_len ) {
			max_len = E - B + 1;
		}
	}

	cur_seq = (char *) ckalloc(sizeof(char) * (max_len+1));
	for( i = 0; i < num_genes; i++ ) {
		if( genes[i].exonCount >= 2 ) {
			strcpy(chr_name, genes[i].sname);
			B = genes[i].txStart;
			E = genes[i].txEnd;
			dir = genes[i].strand;
			strcpy(info, genes[i].gname);	

			k = 0;
			if( dir == '+' ) {
				for (j = B; j <= E; j++) {
					cur_seq[k] = s[j];
					k++;
				}
				cur_seq[k] = '\0';
			}
			else {
				k = 0;
				for (j = E; j >= B; j--) {
					cur_seq[k] = compl[s[j]];
					k++;
				}
				cur_seq[k] = '\0';
			}

			is_correct_splicing = true;
			is_correct_splicing = check_introns(genes, i, exons, cur_seq, k);
			if( is_correct_splicing == false ) {
				if( no_branchpoint == false ) {
					genes[i].type = REDUN;
				}
			}
			else {
				if( no_branchpoint == true ) {
					genes[i].type = REDUN;
				}
			}
		}
	}

  num_genes = rm_redun_genes(genes, 0, num_genes-1);
  write_in_gff(genes, num_genes, exons, num_exons);
	
	free(cur_seq);
	for( i = 0; i < NUM_BP_SEQ; i++ ) free(branchpoints[i]);
	free(branchpoints);
	seq_close(sf);
	return EXIT_SUCCESS;
}

bool check_introns(struct g_list *genes, int gid, struct exons_list *exons, char *cur_seq, int seq_len)
{
	int sid = 0, eid = 0;
	int i = 0, j = 0;
	int b = 0, e = 0;
	int cur_b = 0, cur_e = 0;
	char *temp_seq;
	bool is_correct_splicing = true;
	int num_skips = 0;
	int num_nu = 0;
	int last_true = -1;
	int num_loops = 0, num_init_exons = 0;

	if (seq_len > 0 ) {
		temp_seq = (char *) ckalloc( seq_len * sizeof(char));

		b = genes[gid].txStart;
		e = genes[gid].txEnd;
	
		sid = genes[gid].cdsStart;
		eid = genes[gid].cdsEnd;
	
		num_init_exons = eid - sid + 1;
		i = sid;
		last_true = -1;
		while( (i < eid) && (is_correct_splicing == true) ) {
			if( genes[gid].strand == '+' ) {
				cur_b = exons[i].reg.upper - b + 1;
				cur_e = exons[i+1].reg.lower - b; 
			}
			else {
				if( exons[i].reg.upper < exons[i+1].reg.lower ) {
					cur_b = e - exons[i+1].reg.lower + 1; 
					cur_e = e - exons[i].reg.upper;
				}
				else {
					cur_b = e - exons[i].reg.lower + 1; 
					cur_e = e - exons[i+1].reg.upper; 
				}
			}
			
			for( j = cur_b; j <= cur_e; j++ ) {
				temp_seq[j-cur_b] = cur_seq[j];
			}
			temp_seq[j] = '\0';

			if( check_strstr(temp_seq, branchpoints, NUM_BP_SEQ) == false ) {
				if( i >= (eid-1) ) {
					is_correct_splicing = false;	
				}
				else {
					j = i+1;
					num_skips = 1;	
					num_nu = width(exons[j].reg)+1;
					while( (j <= (eid-1)) && (num_nu % 3 != 0) ) {
						j++;
						num_skips++;
						if( j <= eid ) {	
							num_nu = num_nu + width(exons[j].reg) + 1;	
						}
					}

					if(i < (eid-num_skips)) {
						for( j = (i+1); j <= (eid-num_skips); j++ ) {
							assign_exons(exons, j, exons[j+num_skips]);
						}
						eid = eid - num_skips;
						i--;
						genes[gid].cdsEnd = eid;
					}
					else {
						if( last_true != -1 ) {
							i = last_true;
							j = i+1;
							num_skips = 1;	
							num_nu = width(exons[j].reg)+1;
							while( (j <= (eid-1)) && (num_nu % 3 != 0) ) {
								j++;
								num_skips++;
								if( j <= eid ) {	
									num_nu = num_nu + width(exons[j].reg) + 1;	
								}
							}

							if(i < (eid-num_skips)) {
								for( j = (i+1); j <= (eid-num_skips); j++ ) {
									assign_exons(exons, j, exons[j+num_skips]);
								}
								eid = eid - num_skips;
								genes[gid].cdsEnd = eid;
							}
							else {
								is_correct_splicing = false;	
							}
						}
						else {
							is_correct_splicing = false;	
						}
					}
				}
			}
			else {
				last_true = i;
			}
			i++;
			num_loops++;
			if( num_loops >= num_init_exons ) {
				fatalf("too many loops run: %d\n", num_loops);
			}
		}
		free(temp_seq);
	}
	else {
		fatalf("sequence length of %s: %d\n", genes[gid].gname, seq_len);
	}	

	return(is_correct_splicing);
}

bool check_strstr(char *cur_seq, char **branchpoints, int num_bp_seq)
{
	bool res = false;
	int i = 0;

	for( i = 0; i < num_bp_seq; i++ ) {
		if( strstr(cur_seq, branchpoints[i]) != 0 ) {
			res = true;
		}
	}

	return(res);
}
