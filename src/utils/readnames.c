#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define DEFAULT 0
#define SOAP_UNMAPPED 1
#define SOAP_MAP 2

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	char head[1000];
	char name[1000];
	char prev[1000];
	char mode[100];
	int len = 0;
	int num_reads = 0;
	int flag = DEFAULT;

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(name, "");
	strcpy(prev, "");
	strcpy(mode, "");
	if ( argc == 2 ) {
		flag = DEFAULT; // DEFAULT mode is FASTQ
	}
	else if( argc != 3 ) {
		printf("readnames unmapped_soap_file flag(SOAP_UNMAPPED, SOAP_MAP, or FASTQ)\n");
		printf("Or readnames unmapped_soap_file flag(SOAP_UNMAPPED, SOAP_MAP, or FASTQ) number_of_chars(to ignore at the end of the readname in the input)\n");
		return EXIT_FAILURE;
	}
	else {
		strcpy(mode, argv[2]);
		if( strcmp(mode, "SOAP_MAP") == 0 ) {
			flag = SOAP_MAP;
		}
		else if( strcmp(mode, "SOAP_UNMAPPED") == 0 ) {
			flag = SOAP_UNMAPPED;
		}
		else if( strcmp(mode, "FASTQ") == 0 ) {
			flag = DEFAULT;
		}
		else {
			printf("unsupported flag (SOAP_UNMAPPED, SOAP_MAP, or FASTQ)\n");
			return EXIT_FAILURE;
		}
	}

	if(!(f = fopen(argv[1], "r"))) {
		printf("no file %s exists\n", argv[1]);
		return EXIT_FAILURE;
	}

	while(fgets(buf, 10000, f))
	{
		if( sscanf(buf, "%s %*s", head) != 1 ) {
			printf("bad format in %s\n", buf);
			return EXIT_FAILURE;
		}
		else 
		{
			if(((flag == SOAP_UNMAPPED) && (head[0] == '>')) || ((flag == DEFAULT) && (head[0] == '@'))|| (flag == SOAP_MAP)) {
				len = strlen(head);
				strcpy(name, head);
				if( head[len-2] == '/' ) {
					name[len-2] = '\0';
				}

				if( (num_reads != 0) && (strcmp(name, prev) == 0) ) {}
				else {
					if( flag == SOAP_MAP ) {
						if( head[len-2] == '/' ) {
							printf("%s/1\n", name); 	
							printf("%s/2\n", name); 	
						}
						else {
							printf("%s\n", name); 	
						}
					}
					else {
						if( head[len-2] == '/' ) {
							printf("%s/1\n", name+1); 	
							printf("%s/2\n", name+1); 	
						}
						else {
							printf("%s\n", name+1); 	
						}	
					}
				}
				strcpy(prev, name);
				num_reads++;
			}
			else if( (flag == SOAP_UNMAPPED) || (flag == DEFAULT) ) {}
		}
	}

	
	if(((flag == SOAP_UNMAPPED) && (head[0] == '>')) || ((flag == DEFAULT) && (head[0] == '@'))|| (flag == SOAP_MAP)) {
		if( (num_reads != 0) && (strcmp(name, prev) == 0) ) {}
		else {
			if( head[len-2] == '/' ) {
				printf("%s/1\n", name+1); 	
				printf("%s/2\n", name+1); 	
			}
			else {
				printf("%s\n", name+1); 	
			}
		}
	}
	fclose(f);

	return EXIT_SUCCESS;
}
