#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE 1 
#define FALSE 0

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500], scf_name[500];
	int b, e;
	char mode[10];
	int num_lines = 0;
	int cur_num = 1;

	if( argc != 3 ) {
		printf("ext_loc_info loc_file mode\n");
		return EXIT_FAILURE;
	}

	strcpy(buf, "");
	strcpy(name, "");
	strcpy(scf_name, "");
	strcpy(mode, argv[2]);

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f)) num_lines++;

	fseek(f, 0, SEEK_SET);
	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>') 
		{
			if( sscanf(buf+2, "%d %d %s %s", &b, &e, name, scf_name) == 4 ) {
				if( mode[0] == 'b' ) printf("> %d %d %s %s\n", b-3, e+3, name, scf_name);
				else if( mode[0] == 'M') printf("> %d %d %s %s\n", b-3, e, name, scf_name);
				else if( mode[0] == 'P') printf("> %d %d %s %s\n", b, e+3, name, scf_name);
				else printf("%s", buf);
			}
			else if( sscanf(buf+2, "%d %d %s", &b, &e, name) != 3 )
			{
				printf("bad format in %s\n", buf); 
				return EXIT_FAILURE;
			}
			else
			{
				if( mode[0] == 'b' ) printf("> %d %d %s\n", b-3, e+3, name);
				else if( mode[0] == 'M') printf("> %d %d %s\n", b-3, e, name);
				else if( mode[0] == 'P') printf("> %d %d %s\n", b, e+3, name);
				else printf("%s", buf);
			}
		}
		else if( cur_num == 2 )
		{
			if( sscanf(buf, "%d %d", &b, &e) != 2 )
			{
				printf("bad format in %s\n", buf); 
				return 1;
			}
			else 
			{
				if( mode[0] == 'b' ) printf("%d %d\n", b-3, e);
				else if( mode[0] == 'M') printf("%d %d\n", b-3, e);
				else printf("%s", buf);
			}
		}
		else if( cur_num == num_lines )
		{
			if( sscanf(buf, "%d %d", &b, &e) != 2 )
			{
				printf("bad format in %s\n", buf); 
				return 1;
			}
			else 
			{
				if( mode[0] == 'b' ) printf("%d %d\n", b, e+3);
				else if( mode[0] == 'P') printf("%d %d\n", b, e+3);
				else printf("%s", buf);
			}
		}
		else printf("%s", buf);

		cur_num++;
	}

	return 0;
}
