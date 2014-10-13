#include "main.h"
#include "util.h"
#include "bed.h"
#include "util_i.h"

#define BUF_LEN 10000
#define NUM_REG_COLUMNS 9

void print_rest_columns(char *info1, char *info2, int num_strains);
bool is_indels(char *alt);
int main(int argc, char *argv[])
{
	FILE *f, *g;
	char buf[BUF_LEN];
	char other_buf[BUF_LEN];
	int i = 0;
	char chr1[CHR_LEN], ref1[CHR_LEN], alt1[FEATURE_LEN], info1[MAX_LEN];
	int b1 = 0, b2 = 0;
	char chr2[CHR_LEN], ref2[CHR_LEN], alt2[FEATURE_LEN], info2[MAX_LEN];
	float qual1 = (float)0;
	char **strain_names;
	int num_strains = 0;
	bool is_end = false;
	char *token, *end_token;
	bool is_header = false;
	bool is_ref = false;
	int num_ch = 0;

	strcpy(buf, "");
	strcpy(chr1, "");
	strcpy(ref1, "");
	strcpy(alt1, "");
	strcpy(info1, "");
	strcpy(other_buf, "");
	strcpy(chr2, "");
	strcpy(ref2, "");
	strcpy(alt2, "");
	strcpy(info2, "");

	if( (argc != 4) && (argc != 5) ) {
		printf("vcf2db_snp bcftools_file VarScan_file ref/noref (header)\n");
		return EXIT_FAILURE;
	}
	else {
		if( strcmp(argv[3], "ref") == 0 ) {
			is_ref = true;
		}
		else if( strcmp(argv[3], "noref") == 0 ) {
			is_ref = false;
		}
		else fatalf("illegal argument: %s is not ref/noref\n", argv[3]);

		if( argc == 4 ) {
			is_header = false;
		}
		else if( strcmp(argv[4], "header") == 0 ) {
			is_header = true;
		}
		else {
			fatalf("vcf2db_snp bcftools_file VarScan_file (header)\n");
		}

		if(!(f = fopen(argv[1], "r"))) {
			fatalf("no file %s exists\n", argv[3]);
		}
		else {
			while(fgets(buf, BUF_LEN, f) && (buf[0] == '#') && (buf[1] == '#') ) {}	

			if( (buf[0] == '#') && (buf[1] == '#') ) {
				fatalf("head line is missing including strain names in %s\n", argv[1]);
			}
			strcpy(other_buf, buf);
			token = strtok(buf, " \t");
			while( token != NULL ) {
				token = strtok(NULL, " \t");
				num_strains++;
			}
			num_strains = num_strains - NUM_REG_COLUMNS; // exclude regular column names
			strain_names = (char **) ckalloc(num_strains * sizeof(char *));
			for( i = 0; i < num_strains; i++ ) {
				strain_names[i] = (char *) ckalloc(FEATURE_LEN * sizeof(char));
				strcpy(strain_names[i], "");
			}
			
			i = 0;
			token = strtok_r(other_buf, " \t\n", &end_token);
			while( token != NULL ) {
				if( i >= (num_strains + NUM_REG_COLUMNS )) {
					fatalf("Exceed counted number %d\n", num_strains);
				}

				if( i >= NUM_REG_COLUMNS ) strcpy(strain_names[i-NUM_REG_COLUMNS], token);
				token = strtok_r(NULL, " \t\n", &end_token);
				i++;
			}
		}

		if(!(g = fopen(argv[2], "r"))) {
			printf("no file %s exists\n", argv[2]);
			return EXIT_FAILURE;
		}
	}

	if( is_header == true ) {
		if( is_ref == true ) {
			printf("#{\"column_names\":[\"chr\",\"pos\",\"A\",\"B\",\"qual\",\"ref\",\"rPos\",\"rnuc\"");
			num_ch = num_ch + 64;
		}
		else {
			printf("#{\"column_names\":[\"chr\",\"pos\",\"A\",\"B\",\"qual\"");
			num_ch = num_ch + 43;
		}

		for( i = 0; i < num_strains; i++ ) {
			printf(",\"%dA\",\"%dB\",\"%dG\",\"%dQ\"", i+1, i+1, i+1, i+1);
			num_ch = num_ch + 16 + (4 * ((int)((i+1)/10) + 1));
			if( num_ch >= 190 ) {
				printf("\n#");
				num_ch = 0;
			}
		}
		if( num_ch > 150 ) printf("\n#");
		printf(",\"prim\",\"rflp\"],\"dbkey\":\"?\",\"individuals\":[");
		num_ch = 47;
	
		for( i = 0; i < num_strains; i++ ) {
			if( is_ref == true ) {
				if( i == (num_strains - 1) ) 
					printf("[\"%s\",%d]],", strain_names[i], (i*4)+9);
				else printf("[\"%s\",%d],", strain_names[i], (i*4)+9);
			}
			else {
				if( i == (num_strains - 1) ) 
					printf("[\"%s\",%d]],", strain_names[i], (i*4)+6);
				else printf("[\"%s\",%d],", strain_names[i], (i*4)+6);
			}

			num_ch = num_ch + 7 + strlen(strain_names[i]) + (int)(((i*4)+6)/10);
			if( num_ch >= 190 ) {
				printf("\n#");
				num_ch = 0;
			}
		}
		if( num_ch > 160 ) printf("\n#");
		printf("\"pos\":2,\"rPos\":0,\"ref\":0,\"scaffold\":1,\"species\":\"sacCer\"}\n");
	}

	while( (fgets(buf, BUF_LEN, f)) && (buf[0] == '#') ) {}

	if( buf[0] == '#' ) is_end = true;
	i = 0;
	while(is_end == false)
	{
		if( sscanf(buf, "%s %d %*s %s %s %f %*s", chr1, &b1, ref1, alt1, &qual1) != 5 ) {
			fatalf("bad format in the new bed file 1: %s\n", buf);
		}
		else {
			if( fgets(other_buf, BUF_LEN, g) ) {
				if( (sscanf(other_buf, "%s %*s", chr2) == 1) && (strcmp(chr2, "Chrom") == 0) ) {
					fgets(other_buf, BUF_LEN, g);
				}

				if( sscanf(other_buf, "%s %d %s %s %*s", chr2, &b2, ref2, alt2) != 4 ) {
					fatalf("bad format in the new bed file 2: %s\n", other_buf);
				}
				else {
					if( strcmp(chr1, chr2) != 0 ) {
	;					fatalf("different chromosomes: %s - %s\n", chr1, chr2);
					}	
					else {
						if( b1 != b2 ) {
							while( (is_end == false) && (b1 != b2) ) {
								while( (is_end == false) && (b1 > b2) ) {
									if( fgets(other_buf, BUF_LEN, g) ) {
										if( sscanf(other_buf, "%s %d %s %s %*s", chr2, &b2, ref2, alt2) != 4 ) {
											fatalf("bad format in the new bed file 3: %s\n", other_buf);
										}
									}
									else is_end = true;	
								}

								while( (is_end == false) && (b1 < b2) ) {
									if( fgets(buf, BUF_LEN, f) ) {
										if( sscanf(buf, "%s %d %*s %s %s %f %*s", chr1, &b1, ref1, alt1, &qual1) != 5 ) {
											fatalf("bad format in the new bed file 4: %s\n", buf);
										}
									}
									else is_end = true;	
								}
							}
						}
						else if( (strcmp(chr1, chr2) == 0) && (b1 == b2)  ) {
							if( (is_indels(alt1) == false) && (is_indels(alt2) == false) ) {
								if( is_ref == true ) {
									printf("%s\t%d\t%s\t%s\t%.0f\t%s\t%d\t%s\t", chr1, b1, ref1, alt1, qual1, chr1, b1, ref1); 	
								}
								else {
									printf("%s\t%d\t%s\t%s\t%.0f\t", chr1, b1, ref1, alt1, qual1); 	
								}
								strcpy(info1, buf);
								strcpy(info2, other_buf);
								print_rest_columns(info1, info2, num_strains);
								printf("0\t0\n");
							}
						}
						else if( is_end != true ) {
							fatalf("unexpected run: %s - %s\n", buf, other_buf);
						} 
					}
				}
			}
			else {
				is_end = true;
			}
		}
		
		if( !fgets(buf, BUF_LEN, f) ) {
			is_end = true;
		}
	}

	for( i = 0; i < num_strains; i++ )
	{
		free(strain_names[i]);
	}
	free(strain_names);
	fclose(g);
	fclose(f);

	return EXIT_SUCCESS;
}

void print_rest_columns(char *info1, char *info2, int num_strains)
{
	char *token, *token2;
	char *chr1, *chr2;
	int num_chr1 = 0, num_chr2 = 0, qual = 0;
	int i = 0, j = 0, k = 0;
	struct gd_snp *strains_info;
	char *end_token, *end_field, *end_chr; 
	char temp[FEATURE_LEN], temp2[FEATURE_LEN];
	
	strains_info = (struct gd_snp *) ckalloc(num_strains * sizeof(struct gd_snp));
	token = strtok_r(info1, " \t", &end_token);
	while( (token != NULL) && (i < num_strains) ) {
		if( k < 9 ) {
			token = strtok_r(NULL, " \t", &end_token);
		}
		else {
			strcpy(temp, token);
			token2 = strtok_r(temp, ":", &end_field);
			j = 0;
			while( token2 != NULL ) {
				if( j == 0 ) {
					strcpy(temp2, token2);
					chr1 = strtok_r(temp2, "/", &end_chr);
					num_chr1 = atoi(chr1);
					chr2 = strtok_r(NULL, "/", &end_chr);
					num_chr2 = atoi(chr2);
				}
				else if( j == 2 ) {
					qual = atoi(token2);
				}
				token2 = strtok_r(NULL, ":", &end_field);
				j++;
			}
			strains_info[i].num_chr = num_chr1 + num_chr2;
			strains_info[i].qual = qual;
			token = strtok_r(NULL, " \t", &end_token);
			i++;
		}
		k++;
	}

	i = 0;
	k = 0;
	token = strtok_r(info2, " \t", &end_token);
	while( (token != NULL ) && (i < num_strains) ) {
		if( k < 10 ) {
			token = strtok_r(NULL, " \t", &end_token);
		}
		else {	
			strcpy(temp, token);
			token2 = strtok_r(temp, ":", &end_field);
			j = 0;
			while( (token2 != NULL) && (j < 6) ) 
			{
				if( j == 2 ) {
					strains_info[i].reads1 = atoi(token2);
				}
				else if( j == 3 ) {
				 	strains_info[i].reads2 = atoi(token2);
				}
				token2 = strtok_r(NULL, ":", &end_field);
				j++;
			}
			token = strtok_r(NULL, " \t", &end_token);
			i++;
		}
		k++;
	}

	for( i = 0; i < num_strains; i++ ) {
		printf("%d\t%d\t%d\t%d\t", strains_info[i].reads1, strains_info[i].reads2, strains_info[i].num_chr, strains_info[i].qual);
	}
	
	free(strains_info);
}

bool is_indels(char *alt)
{
	char temp[FEATURE_LEN];
	char *token;
	bool res = false;

	strcpy(temp, alt);
	token = strtok(temp, ",");
	while( (token != NULL) && (res == false) ) {
		if(strlen(token) >= 2) res = true;	
		token = strtok(NULL, ",");
	}
	return(res);	
}
