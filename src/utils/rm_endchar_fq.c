#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int len1 = 0, len2 = 0, len = 0;
	int diff = 0;
	int num = 0;

	strcpy(buf, "");
	if( argc != 2 ) {
		printf("rm_endchar fastq_file\n");
		return EXIT_FAILURE;
	}

	if(!(f = fopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	while(fgets(buf, 10000, f))
	{
		num++;
		if( (num % 2) == 1 ) {
			len = strlen(buf);
			if( buf[len-1] == '\n' ) {
				printf("%s", buf);
			}
			else printf("%s\n", buf);
		}
		else if( (num % 4) == 2 ) {
			len1 = strlen(buf);
			if( buf[len1-1] == '\n' ) {
				printf("%s", buf);
			}
			else printf("%s\n", buf);
		}
		else {
			len2 = strlen(buf);
			diff = 0;
			if( len2 > len1 ) {
				diff = len2 - len1;
			}
			else if( len2 < len1 ) {
				printf("missing quality scores %d-%d\n", len2, len1);
				return EXIT_FAILURE;
			}

			if( buf[len2-1] == '\n' ) {
				buf[len2-diff-1] = '\0';
			}
			else {
				buf[len2-diff] = '\0';
			}
			printf("%s\n", buf);
			
			num = 0;
		}
	}

	return EXIT_SUCCESS;
}
