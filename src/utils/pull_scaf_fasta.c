#include <stdio.h>
#include <string.h>

#define true 1
#define false 0

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500], scaf_name[500], loc[500];
	int is_end = 0;
	int is_gname_given = false;
	char head ='\0';

	strcpy(buf, "");
	strcpy(name, "");
	strcpy(scaf_name, "");
	strcpy(loc, "");

	if( (argc == 4) && (strcmp(argv[3], "GENE_NAME") == 0 ))  {
		is_gname_given = true;
	}
	else if( argc != 3 ) {
		printf("pull_scaf_fasta fasta scaf_name\n");
		return 1;
	}

	f = fopen(argv[1], "r");

	while((is_end == 0) && fgets(buf, 500, f))
	{
		if((buf[0] == '>') || (buf[0] == '<'))
		{
			head = buf[0];
			if( is_gname_given == true ) {
				if( sscanf(buf+1, "%s %s %s %*s", scaf_name, loc, name) != 3 ) {
					printf("bad format in %s\n", buf);
					return 1;
				}
				else {
					if( strcmp(name, argv[2]) == 0 ) is_end = 1;
				}
			}
			else {
				if( sscanf(buf+1, "%s %*s", name) != 1 ) {
					printf("bad format in %s\n", buf);
					return 1;
				}
				else {
					if( strcmp(name, argv[2]) == 0 ) is_end = 1;
				}
			}
		}
	}

	if( is_gname_given == true ) {
		printf("%c%s %s %s\n", head, name, scaf_name, loc);
	}
	else printf("%s", buf);

	while(fgets(buf, 500, f) && (buf[0] != '>')) printf("%s", buf);

	return 0;
}
