#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500];
	char next[500];
	int len = 0;
	int cutoff = 0;
	int scf_size = 0;

	if( argc != 3 ) {
		printf("trim_short_contigs seq_file(fasta) cut_off\n");
		return EXIT_FAILURE;
	}
	else {
		cutoff = atoi(argv[2]);
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>') 
		{
			if( sscanf(buf+1, "%s %d", next, &scf_size) != 2 ) {
				printf("not SGA-scaffolds format %s", buf);
				return EXIT_FAILURE;
			}
			else 
			{
				if( scf_size >= cutoff ) {
					printf(">%s", buf);
				}
			}
		}
		else {
			if( scf_size >= cutoff ) {
				len = strlen(buf);
				if( buf[len-1] == '\n' ) printf("%s", buf);
				else printf("%s\n", buf);
			}
		}
	}

	fclose(f);
	return EXIT_SUCCESS;
}
