// input: NCBI SRA FASTQ format

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define READ_LEN 1000
#define READ1 1
#define READ2 2

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int len = 0, num = 0;
	char head[READ_LEN];
	char name[READ_LEN];
	char part[READ_LEN];
	int type = 0;
	char read_num = '\0';

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(name, "");
	strcpy(part, "");

	if( argc != 2 ) {
		printf("fastq_head fastq_file\n");
		return EXIT_FAILURE;
	}

	if(!(f = fopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	while(fgets(buf, 10000, f))
	{
		num++;
		if( (num % 4) == 1 ) {
			len = strlen(buf);
			if( buf[len-1] == '\n' ) {
				read_num = buf[len-2];
				if( read_num == '1' ) {
					type = READ1;
					strcpy(head, buf);
					head[len-3] = '\0';
				}
				else if( read_num == '2' ) {
					type = READ2;
					strcpy(head, buf);
					head[len-3] = '\0';
				}
				else {
					printf("head line in wrong format: %s\n", buf);
					return EXIT_FAILURE;
				}
			}
			else {
				read_num = buf[len-1];
				if( read_num == '1' ) {
					type = READ1;
					strcpy(head, buf);
					head[len-2] = '\0';
				}
				else if( read_num == '2' ) {
					type = READ2;
					strcpy(head, buf);
					head[len-2] = '\0';
				}
				else {
					printf("head line in wrong format: %s\n", buf);
					return EXIT_FAILURE;
				}
			}

			if( read_num == '1' ) {
				if( sscanf(head, "%[^#]#%s", name, part ) != 2 ) {
					printf("head line in wrong format: %s\n", head);
					return EXIT_FAILURE;
				}
				printf("%s 1:%s\n", name, part); 		
			}
			else if( read_num == '2' ) {
				if( sscanf(head, "%[^#]#%s", name, part ) != 2 ) {
					printf("head line in wrong format: %s\n", head);
					return EXIT_FAILURE;
				}

				printf("%s 2:%s\n", name, part); 		
			}
			else {
				printf("head line in wrong format: %s\n", buf);
				return EXIT_FAILURE;
			}
		}
		else {
			printf("%s", buf);
		}
	}

	fclose(f);
	return EXIT_SUCCESS;
}
