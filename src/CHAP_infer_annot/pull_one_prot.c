#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500];
	int is_end = 0;

	strcpy(buf, "");
	strcpy(name, "");

	if( argc != 3 ) {
		printf("pull_one_proc protein_seq human_gname\n");
		return 1;
	}

	f = fopen(argv[1], "r");

	while((is_end == 0) && fgets(buf, 500, f))
	{
		if(buf[0] == '>')
		{
			if( sscanf(buf+1, "%s %*s", name) != 1 ) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else {
				if( strcmp(name, argv[2]) == 0 ) is_end = 1;
			}
		}
	}

	printf("%s", buf);
	while(fgets(buf, 500, f) && (buf[0] != '>')) printf("%s", buf);

	return 0;
}
