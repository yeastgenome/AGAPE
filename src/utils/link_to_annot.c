#include "main.h"
#include "util.h"
#include "util_i.h"
#include "util_annot.h"
#include "seq.h"

int b[1000], e[1000], N;
char compl[128];

char GeneticCode[4][4][4] = {'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'};

int find_overlap_gene(char *chr_name, int b, struct exons_list *exons, int num_cds)
{
	int id = -1;
	int i = 0;
	bool is_in = false;
	
	while( (i < num_cds) && (is_in == false) ) {
		if( (strcmp(chr_name, exons[i].chr) == 0 ) && (in(b, exons[i].reg) == true) ) 
		{
			id = i;
			is_in = true;
		}
		i++;
	}

	return(id);
}

char dna2oneaa(char *dna)
{
	int index = 0;
	int frame1 = 0, frame2 = 0, frame3 = 0;

	if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u') frame1 = 0;
	else if(dna[index] == 'C' || dna[index] == 'c')
		frame1 = 1;
	else if(dna[index] == 'A' || dna[index] == 'a')
		frame1 = 2;
    else if(dna[index] == 'G' || dna[index] == 'g')
      frame1 = 3;
    else {
      printf("\nWrong DNA format\n");
      return EXIT_FAILURE;
    }

    index++;

    if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u')
      frame2 = 0;
    else if(dna[index] == 'C' || dna[index] == 'c')
      frame2 = 1;
    else if(dna[index] == 'A' || dna[index] == 'a')
      frame2 = 2;
    else if(dna[index] == 'G' || dna[index] == 'g')
      frame2 = 3;
    else {
      printf("\nWrong DNA format\n");
      return EXIT_FAILURE;
    }

    index++;

    if(dna[index] == 'T' || dna[index] == 't' || dna[index] == 'U' || dna[index] == 'u')
      frame3 = 0;
    else if(dna[index] == 'C' || dna[index] == 'c')
      frame3 = 1;
    else if(dna[index] == 'A' || dna[index] == 'a')
      frame3 = 2;
    else if(dna[index] == 'G' || dna[index] == 'g')
      frame3 = 3;
    else {
      printf("\nWrong DNA format\n");
      return EXIT_FAILURE;
    }

		return(GeneticCode[frame1][frame2][frame3]);
}

int main(int argc, char *argv[])
{
	SEQ *sf;
	uchar *s;
	FILE *f;
	char buf[10000];
	char head[MAX_LEN];
	char cur[LEN_NAME], chr_name[LEN_NAME], annot[LEN_NAME], gname[LEN_NAME], filter[LEN_NAME];
	int gid = -1;
	int rid = -1;
	int i = 0;
	int b = 0, e = 1, num_cds = 0;
	char dir[3];
	struct exons_list *exons;
	char annot_name[LEN_NAME];
	float qual = (float)0;
	char ref[LEN_NAME], alt[LEN_NAME];
	int rest = 0;
	char codon[4], alt_codon[4];
	char aa1 = '\0', aa2 = '\0';
	int num_rmsk = 0;
	struct exons_list *rmsk;
	int num_snps = 0, num_pass = 0, num_filter = 0, num_coding1 = 0, num_syn1 = 0, num_non1 = 0, num_repeats1 = 0, num_coding_repeats1 = 0;
	int num_coding = 0, num_syn = 0, num_non = 0, num_repeats = 0, num_coding_repeats = 0;
	bool is_num_print = false;

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(cur, "");
	strcpy(chr_name, "");
	strcpy(annot, "");
	strcpy(gname, "");
	strcpy(annot_name, "");
	strcpy(ref, "");
	strcpy(alt, "");
	strcpy(codon, "");
	strcpy(alt_codon, "");
	strcpy(dir, "");
	codon[3] = '\0';
	alt_codon[3] = '\0';

	if( argc != 7 ) {
		printf("link_to_annot vcf_file gff_file seq_file annot_type(exon, gene, ...) rmsk_file print_mode(NUM or SITES)\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = ckopen(argv[2], "r"))) {
			printf("no file %s exists\n", argv[2]);
			return EXIT_FAILURE;
		}

		strcpy(annot_name, argv[4]);
		if( strcmp(annot_name, "exon") != 0 ) {
			fatalf("seq file is required only when the annot type is exon, but %s here\n", annot_name);
		}
		sf = seq_get(argv[3]);
		s = SEQ_CHARS(sf) - 1;
		if( strcmp(argv[6], "NUM") == 0 ) {
			is_num_print = true;
		}
		else if( strcmp(argv[6], "SITES") == 0 ) {
			is_num_print = false;
		}
		else {
			fatalf("unsupported print option: %s\n", argv[6]);
		}
	}

  compl['a'] = compl['A'] = 'T';
  compl['c'] = compl['C'] = 'G';
  compl['g'] = compl['G'] = 'C';
  compl['t'] = compl['T'] = 'A';

	while(fgets(buf, 10000, f))
	{
		if( (buf[0] == '#') || (buf[0] == '>') ) {}
		else if( sscanf(buf, "%*s %*s %s %d %d %*s", annot, &b, &e) != 3 ) {
			fatalf("line in wrong gff format: %s\n", buf);
		}
		else {
			if( strcmp(annot, annot_name) == 0 ) {
				num_cds++;
			}
		}
	}

	if( num_cds > 0 ) exons = (struct exons_list *) ckalloc(num_cds * sizeof(struct exons_list));

	initialize_exons_list(exons, 0, num_cds);

	fseek(f, 0, SEEK_SET);

	i = 0;
	
	while(fgets(buf, 10000, f))
	{
		if( (buf[0] == '#') || (buf[0] == '>') ) {}
		else if( sscanf(buf, "%s %*s %s %d %d %*s %s %*s %s", chr_name, annot, &b, &e, dir, cur) != 6 ) {
			fatalf("line in wrong gff format: %s\n", buf);
		}
		else {
			if( strcmp(annot, annot_name) == 0 ) {
				get_gene_name(cur, gname);
				strcpy(exons[i].name, gname);
				exons[i].reg = assign_I(b, e);
				exons[i].dir = dir[0];
				strcpy(exons[i].chr, chr_name);
				i++;
			}	
		}
	}

	if( i != num_cds ) {
		fatalf("%s counting error: %d - %d\n", annot_name, num_cds, i);
	}
	fclose(f);

	if(!(f = ckopen(argv[5], "r"))) {
		fatalf("%s file not found\n", argv[5]);
	}

	rmsk = 0;
	while(fgets(buf, 10000, f))
	{
		if( (buf[0] == '#') || (buf[0] == '>') ) {}
		else if( sscanf(buf, "%*s %*s %s %d %d %*s", annot, &b, &e) != 3 ) {
			fatalf("line in wrong gff format: %s\n", buf);
		}
		else {
			num_rmsk++;
		}
	}

	if( num_rmsk > 0 ) rmsk = (struct exons_list *) ckalloc(num_rmsk * sizeof(struct exons_list));

	initialize_exons_list(rmsk, 0, num_rmsk);

	fseek(f, 0, SEEK_SET);

	i = 0;
	
	while(fgets(buf, 10000, f))
	{
		if( (buf[0] == '#') || (buf[0] == '>') ) {}
		else if( sscanf(buf, "%s %*s %s %d %d %*s %s %*s %s", chr_name, annot, &b, &e, dir, cur) != 6 ) {
			fatalf("line in wrong gff format: %s\n", buf);
		}
		else {
			strcpy(rmsk[i].name, annot);
			rmsk[i].reg = assign_I(b, e);
			rmsk[i].dir = dir[0];
			strcpy(rmsk[i].chr, chr_name);
			i++;
		}
	}

	if( i != num_rmsk ) {
		fatalf("%s counting error: %d - %d\n", annot_name, num_cds, i);
	}
	fclose(f);

	if(!(f = ckopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	i = 0;
	while(fgets(buf, 10000, f))
	{
		if( buf[0] != '#' ) {
			num_snps++;
			if( sscanf(buf, "%s %d %*s %s %s %f %s %*s", chr_name, &b, ref, alt, &qual, filter) != 6 ) {
				fatalf("bad format in %s\n", buf);
			}
			else 
			{
				if( strstr(filter, "PASS") == 0 ) {
					num_pass++;
				}
				else if( strstr(filter, "filter") == 0 ) {
					num_filter++;
				}

				rid = -1;
				rid = find_overlap_gene(chr_name, b, rmsk, num_rmsk);
				if( rid != -1 ) {
					num_repeats++;
					if( strstr(filter, "filter") == 0 ) {}
					else if( strstr(filter, "PASS") == 0 ) {
						num_repeats1++;
					}
					else {
						fatalf("unexpected filter option: %s\n", filter);
					}
				}

				if( (gid = find_overlap_gene(chr_name, b, exons, num_cds)) != -1 ) {
					num_coding++;
					if( strstr(filter, "PASS") == 0 ) {
						num_coding1++;
					}

					if( ref[0] != s[b] ) {
						fatalf("nucleotides not match: %c - %c\n", alt, s[b]);
					}

					if( exons[gid].dir == '+' ) {
						rest = (b - exons[gid].reg.lower)%3;
						if( rest == 0 ) {
							sprintf(codon, "%c%c%c", s[b], s[b+1], s[b+2]);
							sprintf(alt_codon, "%c%c%c", alt[0], s[b+1], s[b+2]);
						}
						else if( rest == 1 ) {
							sprintf(codon, "%c%c%c", s[b-1], s[b], s[b+1]);
							sprintf(alt_codon, "%c%c%c", s[b-1], alt[0], s[b+1]);
						}
						else {
							sprintf(codon, "%c%c%c", s[b-2], s[b-1], s[b]);
							sprintf(alt_codon, "%c%c%c", s[b-2], s[b-1], alt[0]);
						}
					}
					else if( exons[gid].dir == '-' ) {
						rest = (b - exons[gid].reg.upper)%3;
						if( rest == 0 ) {
							sprintf(codon, "%c%c%c", compl[s[b]], compl[s[b-1]], compl[s[b-2]]);
							sprintf(alt_codon, "%c%c%c", compl[alt[0]], compl[s[b-1]], compl[s[b-2]]);
						}
						else if( rest == 1 ) {
							sprintf(codon, "%c%c%c", compl[s[b+1]], compl[s[b]], compl[s[b-1]]);
							sprintf(alt_codon, "%c%c%c", compl[s[b+1]], compl[alt[0]], compl[s[b-1]]);
						}
						else {
							sprintf(codon, "%c%c%c", compl[s[b+2]], compl[s[b+1]], compl[s[b]]);
							sprintf(alt_codon, "%c%c%c", compl[s[b+2]], compl[s[b+1]], compl[alt[0]]);
						}
					}
					else {
						fatalf("%c unsupported\n", exons[gid].dir);
					}
					aa1 = dna2oneaa(codon);
					aa2 = dna2oneaa(alt_codon);
					
					if( aa1 == aa2 ) {
						num_syn++;
						if( strstr(filter, "filter") == 0) {
						}
						else if( strstr(filter, "PASS") == 0 ) {
							num_syn1++;
						}
						else {
							fatalf("unexpected filter option: %s\n", filter);
						}
					}
					else {
						num_non++;
						if( strstr(filter, "filter") == 0) {
						}
						else if( strstr(filter, "PASS") == 0 ) {
							num_non1++;
						}
						else {
							fatalf("unexpected filter option: %s\n", filter);
						}
					}
					
					if( rid != -1 ) {
						num_coding_repeats++;
						if( strstr(filter, "PASS") == 0 ) {
							num_coding_repeats1++;
						}
					}

					if( is_num_print == false ) {
						if( rid == -1 ) {
							printf("%s\t%d\t%s\t%s\t%f\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t.\n", chr_name, b, ref, alt, qual, filter, exons[gid].name, exons[gid].reg.lower, exons[gid].reg.upper, exons[gid].dir, aa1, aa2);
						}
						else {
							printf("%s\t%d\t%s\t%s\t%f\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t%s\n", chr_name, b, ref, alt, qual, filter, exons[gid].name, exons[gid].reg.lower, exons[gid].reg.upper, exons[gid].dir, aa1, aa2, rmsk[rid].name);
						}
					}
				}
				else {
				}
			}
		}
	}
	
	if( is_num_print == true ) {
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", chr_name, num_snps, num_pass, num_filter, num_coding, num_coding1, num_non, num_syn, num_non1, num_syn1, num_repeats, num_repeats1, num_coding_repeats, num_coding_repeats1);
	}

	if( num_cds > 0 ) {
		free(exons);
	}
	fclose(f);

	return EXIT_SUCCESS;
}
