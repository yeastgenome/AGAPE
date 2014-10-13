#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define TYPE1 1
#define TYPE2 2

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int len = 0, num = 0;
	char end_chars[100];
	char read_num[100];
	char head1[1000], head2[1000];
	int type = 0;

	strcpy(buf, "");
	strcpy(read_num, "");
	strcpy(head1, "");
	strcpy(head2, "");
	strcpy(end_chars, "");

	if( argc != 3 ) {
		printf("edit_fq_head fastq_file read_num\n");
		return EXIT_FAILURE;
	}
	strcpy(read_num, argv[2]);

	if(!(f = fopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	while(fgets(buf, 10000, f))
	{
		num++;
		if( (num % 4) == 1 ) {
			if( sscanf(buf, "%s %s %*s", head1, head2) == 2 ) {
				type = TYPE2;
			}
			else if( sscanf(buf, "%s %*s", head1) != 1 ) {
				type = TYPE1;
				printf("bad format in head %s\n", buf);
				return EXIT_FAILURE;
			}

/*
			if( sscanf(buf, "%s %*s", head1) != 1 ) {
				printf("bad format in head %s\n", buf);
				return EXIT_FAILURE;
			}
*/
			if( type == TYPE2 ) {
				printf("%s#%s/%s\n", head1, head2+2, read_num);
			}
			else if( type == TYPE1 ) {
				printf("%s/%s\n", head1, read_num);
			}
		}
		else {
			len = strlen(buf);
			if( buf[len-1] == '\n' ) {
				printf("%s", buf);
			}
			else printf("%s\n", buf);
		}
	}

	return EXIT_SUCCESS;
}
