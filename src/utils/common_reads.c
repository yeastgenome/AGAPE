// input: fastq1 (big in SRA format) fastq2 (short ending with /)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"

#define READ_LEN 1000
#define READ1 1
#define READ2 2

#define TRUE 0
#define FALSE 1

int is_common_reads(char *name, char **heads, int num_reads);
int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int len = 0, num = 0;
	char head[READ_LEN];
	char name[READ_LEN];
	char part[READ_LEN];
	char **heads2;
	int num_reads2 = 0;
	int i = 0;
	int is_common = FALSE;

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(name, "");
	strcpy(part, "");

	if( argc != 3 ) {
		printf("common_reads fastq_file1 fastq_file2\n");
		return EXIT_FAILURE;
	}

	if(!(f = ckopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	num = 0;
	while( fgets(buf, 10000, f) ) {
		num++;
		if( (num % 4) == 1 ) {
			num_reads2++;
		}
	}
	
	if( num_reads2 > 0 ) {
		heads2 = (char **) ckalloc(sizeof(char *) * num_reads2);
		for( i = 0; i < num_reads2; i++ ) {
			heads2[i] = (char *) ckalloc(sizeof(char) * READ_LEN);
			strcpy(heads2[i], "");
		}
	}

	fseek(f, 0, SEEK_SET);
	num = 0;
	i = 0;
	while( fgets(buf, 10000, f) ) {
		num++;
		if( (num % 4) == 1 ) {
			len = strlen(buf);
			strcpy(heads2[i], buf);
			if( buf[len-1] == '\n' ) {
				heads2[i][len-3] = '\0';
			}
			else {
				heads2[i][len-2] = '\0';
			}
			i++;
		}
	}
	fclose(f);

	if(!(f = ckopen(argv[2], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	num = 0;
	is_common = FALSE;
	while(fgets(buf, 10000, f))
	{
		num++;
		if( (num % 4) == 1 ) {
			if( sscanf(buf, "%s %s", name, part ) != 2 ) {
				printf("head line in wrong format: %s\n", buf);
				return EXIT_FAILURE;
			}
			is_common = is_common_reads(name, heads2, num_reads2);
			if( is_common == TRUE ) {
				printf("%s", buf);
			}
		}
		else {
			if( is_common == TRUE ) {
				printf("%s", buf);
			}
		}
	}

	if( num_reads2 > 0 ) {
		for( i = 0; i < num_reads2; i++ ) {
			free(heads2[i]);
		}
		free(heads2);
	}

	fclose(f);
	return EXIT_SUCCESS;
}

int is_common_reads(char *name, char **heads, int num_reads)
{
	int i = 0;
	int res = FALSE;

	while( (res == FALSE) && (i < num_reads) ) {
		if( strcmp(name, heads[i]) == 0 ) {
			res = TRUE;
		}
		i++;
	}

	return(res);
}
