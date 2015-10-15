#include "main.h"
#include "util.h"
#include "bed.h"
#include "util_i.h"

#define BUF_LEN 10000

int main(int argc, char *argv[])
{
	FILE *f, *g;
	char buf[BUF_LEN];
	struct bed *snps;
	int num_rows = 0;
	int i = 0, j = 0;
	char chr[CHR_LEN], ref[CHR_LEN], alt[CHR_LEN], pass[FEATURE_LEN], info[MAX_LEN];
	int b = 0, e = 1;
	float qual = (float)0;
	char strain_name[FEATURE_LEN];

	strcpy(buf, "");
	strcpy(chr, "");
	strcpy(ref, "");
	strcpy(alt, "");
	strcpy(pass, "");
	strcpy(info, "");
	if( (argc != 3) && (argc != 4) ) {
		printf("merge_bed file1 file2 (or strain_name)\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = fopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}

		if(!(g = fopen(argv[2], "r"))) {
			printf("no file %s exists\n", argv[2]);
			return EXIT_FAILURE;
		}

		if( argc == 4 ) {
			strcpy(strain_name, argv[3]);
		}
	}

	while(fgets(buf, BUF_LEN, g))
	{
		num_rows++;
	}

	fseek(g, 0, SEEK_SET);

	i = 0;
	snps = (struct bed *) ckalloc(num_rows * sizeof(struct bed));
	init_bed(snps, num_rows);
	while(fgets(buf, BUF_LEN, g))
	{
		if( sscanf(buf, "%s %d %d %*s %f %s %s %s %s", chr, &b, &e, &qual, ref, alt, pass, info) != 8 ) {
			fatalf("bad format in the new bed file: %s in %s\n", buf, strain_name);
		}
		else {
			strcpy(snps[i].chr, chr);
			snps[i].reg = assign_I(b, e);
			snps[i].qual = qual;
			strcpy(snps[i].ref, ref);
			strcpy(snps[i].alt, alt);
			strcpy(snps[i].pass, pass);
			strcpy(snps[i].info, info);
			i++;
		}
	}

	if( i != num_rows ) {
		fatalf("counting error: %d - %d\n", i, num_rows);
	}

	j = 0;
	while( fgets(buf, BUF_LEN, f) ) {	
		if( sscanf(buf, "%s %d %d %*s %f %s %s %s %s", chr, &b, &e, &qual, ref, alt, pass, info) != 8 ) {
			fatalf("bad format in the merged file: %s in %s\n", buf, strain_name);
		}
		else 
		{
			if( (j < num_rows) && (strcmp(chr, snps[j].chr) != 0) ) {
				fatalf("different chromosomes: %s - %s\n", chr, snps[j].chr);
			}
			else {
				if( (j >= num_rows) || ((j < num_rows) && (b < snps[j].reg.lower)) ) {
					printf("%s\t%d\t%d\t.\t%.2f\t%s\t%s\t%s\t%s\n", chr, b, e, qual, ref, alt, pass, info);
				}
				else if( b == snps[j].reg.lower ) {
					while( (j < num_rows) && (b == snps[j].reg.lower) ) {
						if( qual < snps[j].qual ) {
							e = snps[j].reg.upper;
							qual = snps[j].qual;
							if( strcmp(ref, snps[j].ref) != 0 ) {
								fatalf("different reference nucleotides: %s - %s\n", ref, snps[j].ref);
							}
							strcpy(alt, snps[j].alt);
							strcpy(pass, snps[j].pass);
							strcpy(info, snps[j].info);
						}
						j++;
					}
					printf("%s\t%d\t%d\t.\t%.2f\t%s\t%s\t%s\t%s\n", chr, b, e, qual, ref, alt, pass, info);
				}
				else { // b > snps[j].reg.lower
					while( (j < num_rows) && (b > snps[j].reg.lower) ) {
						printf("%s\t%d\t%d\t.\t%.2f\t%s\t%s\t%s\t%s\n", snps[j].chr, snps[j].reg.lower, snps[j].reg.upper, snps[j].qual, snps[j].ref, snps[j].alt, snps[j].pass, snps[j].info);
						j++;
					}

					while( (j < num_rows) && (b == snps[j].reg.lower) ) {
						if( qual < snps[j].qual ) {
							e = snps[j].reg.upper;
							qual = snps[j].qual;
							if( strcmp(ref, snps[j].ref) != 0 ) {
								fatalf("different reference nucleotides: %s - %s\n", ref, snps[j].ref);
							}
							strcpy(alt, snps[j].alt);
							strcpy(pass, snps[j].pass);
							strcpy(info, snps[j].info);
						}
						j++;
					}
					printf("%s\t%d\t%d\t.\t%.2f\t%s\t%s\t%s\t%s\n", chr, b, e, qual, ref, alt, pass, info);
				}
			}
		}
	}

	free(snps);
	fclose(g);
	fclose(f);

	return EXIT_SUCCESS;
}
