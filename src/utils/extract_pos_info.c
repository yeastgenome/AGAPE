#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	char head[10];
	char sign[10];
	char scf_name[100];
	int len = 0, pos = 0;

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(sign, "");
	strcpy(scf_name, "");
	if( argc != 2 ) {
		printf("extract_pos_info soap_file\n");
		return EXIT_FAILURE;
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 10000, f))
	{
		if( sscanf(buf, "%*s %*s %*s %*s %s %d %s %s %d %*s", head, &len, sign, scf_name, &pos) != 5 ) {
			printf("bad format in %s\n", buf);
			return EXIT_FAILURE;
		}
		else 
		{
			if(strcmp(head, "a") == 0) {
				printf("1\t%d\t%s\t%s\t%d\n", len, sign, scf_name, pos); 	
			}
			else if(strcmp(head, "b") == 0) {
				printf("2\t%d\t%s\t%s\t%d\n", len, sign, scf_name, pos); 	
			}
		}
	}

	fclose(f);

	return EXIT_SUCCESS;
}
