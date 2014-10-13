#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TRUE 1
#define FALSE 0
#define TH 30

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500], chr[500];
	int b = 0, e = 1;

	strcpy(buf, "");
	strcpy(name, "");
	strcpy(chr, "");

	if( argc != 2 ) {
		printf("get_gene_bound loc_file\n");
		return 1;
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if((buf[0] == '>') || (buf[0] == '<')) {
			sscanf(buf+2, "%s %d %d %s %*s", chr, &b, &e, name);
			if( buf[0] == '>' ) printf("%c %s %d %d %s\n", buf[0], chr, b, e, name);
			else if( buf[0] == '<' ) printf("%c %s %d %d %s (complement)\n", buf[0], chr, b, e, name);
			printf("%d %d\n", b, e);
		}
	}

	fclose(f);

	return EXIT_SUCCESS;
}
