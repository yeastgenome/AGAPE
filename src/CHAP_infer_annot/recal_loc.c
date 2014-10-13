#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500], rest[500];
	int b = 0, e = 1;
	int base = 0;

	if( argc != 3 ) {
		printf("recal_loc loc_file beg_loc\n");
		return 1;
	}

	f = fopen(argv[1], "r");
	base = atoi(argv[2]);

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>') 
		{
			if( sscanf(buf+2, "%d %d %s", &b, &e, name) != 3 ) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else 
			{
				b = b - base;
				e = e - base;
				printf("%c %d %d %s\n", buf[0], b, e, name);
			}
		}
		else if(buf[0] == '<')
		{
			if( sscanf(buf+2, "%d %d %s %s", &b, &e, name, rest) != 4 ) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else 
			{
				b = b - base;
				e = e - base;
				printf("%c %d %d %s %s\n", buf[0], b, e, name, rest);
			}
		}
		else {
			if( sscanf(buf, "%d %d", &b, &e) != 2 ) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else 
			{
				b = b - base;
				e = e - base;
				printf("%d %d\n", b, e);
			}
		}
	}

	fclose(f);

	return 0;
}
