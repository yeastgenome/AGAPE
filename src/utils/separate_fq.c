// input: NCBI SRA FASTQ format

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define VERSION1 1
#define VERSION2 2

#define READ_LEN 1000
#define READ1 1
#define READ2 2

#define PAIR_FQ1 1
#define PAIR_FQ2 2
#define SINGLE_FQ1 3
#define SINGLE_FQ2 4

int main(int argc, char *argv[])
{
	FILE *f;
	FILE *g1, *g2, *g3, *g4;
	char buf[10000];
	int len = 0, num = 0;
	char end_chars[100];
	char head1[READ_LEN], head2[READ_LEN];
	char name1[READ_LEN], name2[READ_LEN], prev_name[READ_LEN];
	char part1[READ_LEN], part2[READ_LEN], prev_part[READ_LEN];
	char line1[READ_LEN], line2[READ_LEN], line3[READ_LEN];
	int type = 0;
	int file_type = 0;
	char read_num = '\0', prev_read_num = '\0';
	int version = 0;

	strcpy(buf, "");
	strcpy(head1, "");
	strcpy(head2, "");
	strcpy(name1, "");
	strcpy(name2, "");
	strcpy(prev_name, "");
	strcpy(prev_part, "");
	strcpy(part1, "");
	strcpy(part2, "");
	strcpy(line1, "");
	strcpy(line2, "");
	strcpy(line3, "");
	strcpy(end_chars, "");

	version = VERSION2;
	if( argc == 7 ) {
		if( strcmp(argv[6], "VERSION1") == 0 ) {
			version = VERSION1;
		}
	}
	else if( argc != 6 ) {
		printf("separate_fq fastq_file file1 file2 file3 file4\n");
		return EXIT_FAILURE;
	}

	if(!(f = fopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	g1 = fopen(argv[2], "w");
	g2 = fopen(argv[3], "w");
	g3 = fopen(argv[4], "w");
	g4 = fopen(argv[5], "w");

	while(fgets(buf, 10000, f))
	{
		num++;
		if( (num % 4) == 1 ) {
			len = strlen(buf);
			
			prev_read_num = read_num;
			strcpy(prev_name, name1);
			strcpy(prev_part, part1);
			if( buf[len-1] == '\n' ) {
				read_num = buf[len-2];
				if( read_num == '1' ) {
					type = READ1;
					strcpy(head1, buf);
					head1[len-3] = '\0';
				}
				else if( read_num == '2' ) {
					type = READ2;
					strcpy(head2, buf);
					head2[len-3] = '\0';
				}
				else {
					printf("head line in wrong format: %s\n", buf);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}
			else {
				read_num = buf[len-1];
				if( read_num == '1' ) {
					type = READ1;
					strcpy(head1, buf);
					head1[len-2] = '\0';
				}
				else if( read_num == '2' ) {
					type = READ2;
					strcpy(head2, buf);
					head2[len-2] = '\0';
				}
				else {
					printf("head line in wrong format: %s\n", buf);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}

			if( read_num == '1' ) {
				if( version == VERSION1 ) {
					strcpy(name1, head1);
				}
				else if( sscanf(head1, "%[^#]#%s", name1, part1 ) != 2 ) {
					printf("head line in wrong format: %s\n", head1);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}
			else if( read_num == '2' ) {
				if( version == VERSION1 ) {
					strcpy(name2, head2);
				}
				else if( sscanf(head2, "%[^#]#%s", name2, part2 ) != 2 ) {
					printf("head line in wrong format: %s\n", head2);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}

			}
			else {
				printf("head line in wrong format: %s\n", buf);
				fclose(g1);
				fclose(g2);
				fclose(g3);
				fclose(g4);
				return EXIT_FAILURE;
			}

			if( read_num == '2') {
				if( (prev_read_num == '1') && (strcmp(name2, name1) == 0) ) {
					if( version == VERSION1 ) {
						fprintf(g1, "%s/1\n", name1); 		
					}
					else {
						fprintf(g1, "%s 1:%s\n", name1, part1); 		
					}
					fprintf(g1, "%s", line1); 		
					fprintf(g1, "%s", line2); 		
					fprintf(g1, "%s", line3); 		
					if( version == VERSION1 ) {
						fprintf(g2, "%s/2\n", name2); 		
					}
					else {
						fprintf(g2, "%s 2:%s\n", name2, part2); 		
					}
					file_type = PAIR_FQ2;
				}
				else if( prev_read_num == '2' ) {
					if( version == VERSION1 ) {
						fprintf(g4, "%s/2\n", name2); 		
					}
					else fprintf(g4, "%s 2:%s\n", name2, part2); 		
					file_type = SINGLE_FQ2;
				}
				else {
					if( version == VERSION1 ) {
						fprintf(g3, "%s/1\n", name1); 		
					}
					else fprintf(g3, "%s 1:%s\n", name1, part1); 		
					fprintf(g3, "%s", line1); 		
					fprintf(g3, "%s", line2); 		
					fprintf(g3, "%s", line3); 		
					if( version == VERSION1 ) {
						fprintf(g4, "%s/2\n", name2); 		
					}
					else fprintf(g4, "%s 2:%s\n", name2, part2); 		
					file_type = SINGLE_FQ2;
				}
			}
			else {
				if( prev_read_num == '1' ) {
					if( version == VERSION1 ) {
						fprintf(g3, "%s/1\n", prev_name); 		
					}
					else fprintf(g3, "%s 1:%s\n", prev_name, prev_part); 		
					fprintf(g3, "%s", line1); 		
					fprintf(g3, "%s", line2); 		
					fprintf(g3, "%s", line3); 		
				}
			}
		}
		else if( (num % 4) == 2 ) {
			if( type == READ1 ) {
				strcpy(line1, buf);
			}
			else if( type == READ2 ) {
				if( file_type == PAIR_FQ2 ) {
					fprintf(g2, "%s", buf);
				}	
				else if( file_type == SINGLE_FQ2 ) {
					fprintf(g4, "%s", buf);
				}
				else {
					printf("format error 2 : %s\n", buf);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}
			else {
				printf("format error 2 : %s\n", buf);
				fclose(g1);
				fclose(g2);
				fclose(g3);
				fclose(g4);
				return EXIT_FAILURE;
			}
		}
		else if( (num % 4) == 3 ) {
			if( type == READ1 ) {
				strcpy(line2, buf);
			}
			else if( type == READ2 ) {
				if( file_type == PAIR_FQ2 ) {
					fprintf(g2, "%s", buf);
				}	
				else if( file_type == SINGLE_FQ2 ) {
					fprintf(g4, "%s", buf);
				}
				else {
					printf("format error 3 : %s\n", buf);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}
			else {
				printf("format error 3 : %s\n", buf);
				fclose(g1);
				fclose(g2);
				fclose(g3);
				fclose(g4);
				return EXIT_FAILURE;
			}
		}
		else {
			if( type == READ1 ) {
				strcpy(line3, buf);
			}
			else if( type == READ2 ) {
				if( file_type == PAIR_FQ2 ) {
					fprintf(g2, "%s", buf);
				}	
				else if( file_type == SINGLE_FQ2 ) {
					fprintf(g4, "%s", buf);
				}
				else {
					printf("format error 4 : %s\n", buf);
					fclose(g1);
					fclose(g2);
					fclose(g3);
					fclose(g4);
					return EXIT_FAILURE;
				}
			}
			else {
				printf("format error 4 : %s\n", buf);
				fclose(g1);
				fclose(g2);
				fclose(g3);
				fclose(g4);
				return EXIT_FAILURE;
			}
		}
	}

	prev_read_num = read_num;
	if( prev_read_num == '1' ) {
		strcpy(prev_name, name1);
		strcpy(prev_part, part1);
		if( version == VERSION1 ) {
			fprintf(g3, "%s/1\n", prev_name); 		
		}
		else fprintf(g3, "%s 1:%s\n", prev_name, prev_part); 		
		fprintf(g3, "%s", line1); 		
		fprintf(g3, "%s", line2); 		
		fprintf(g3, "%s", line3); 		
	}

	fclose(g1);
	fclose(g2);
	fclose(g3);
	fclose(g4);
	fclose(f);
	return EXIT_SUCCESS;
}
