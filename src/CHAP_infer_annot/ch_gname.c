#include <stdio.h>
#include "util.h"

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], scf_name[500];
	int b = 0, e = 1;
	int f_b = 0, f_e = 1;
	int num_loc = 0;
	char cur_dir = '>';

	if( argc == 4 ) {
		strcpy(scf_name, argv[3]);
	}
	else if( argc != 3 ) {
		fatal("ch_gname filename gname\n");
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>')
		{
			if( sscanf(buf+2, "%d %d %*s", &b, &e) != 2 ) {
				fatalf("bad format in %s\n", buf);
			}
			else {
				if( num_loc == 0 ) {
					f_b = b;
					f_e = e;
					cur_dir = '>';
					num_loc++;
				}
				else {
					if( buf[0] == cur_dir )
					{
						if( f_b > b ) f_b = b;
						if( f_e < e ) f_e = e;
					}
					else {
						fatalf("bad format in %s\n", buf);
					}	
				}
			}
		}
		else if( buf[0] == '<')
		{
			if( sscanf(buf+2, "%d %d %*s", &b, &e) != 2 ) {
				fatalf("bad format in %s\n", buf);
			}
			else {
				if( num_loc == 0 ) {
					f_b = b;
					f_e = e;
					cur_dir = '<';
					num_loc++;
				}
				else {
					if( buf[0] == cur_dir )
					{
						if( f_b > b ) f_b = b;
						if( f_e < e ) f_e = e;
					}
					else {
						fatalf("bad format in %s\n", buf);
					}
				}
			}
		}
	}

	fseek(f, 0, SEEK_SET);

	if( argc == 3 ) {
		if( cur_dir == '>' ) printf("> %d %d %s\n", f_b, f_e, argv[2]);
		else if( cur_dir == '<' ) printf("< %d %d %s (complement)\n", f_b, f_e, argv[2]);
	}
	else {
		if( cur_dir == '>' ) printf("> %d %d %s %s\n", f_b, f_e, argv[2], scf_name);
		else if( cur_dir == '<' ) printf("< %d %d %s %s (complement)\n", f_b, f_e, argv[2], scf_name);
	}

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>')
		{
			if( sscanf(buf+2, "%d %d %*s", &b, &e) != 2 ) {
				fatalf("bad format in %s\n", buf);
			}
		}
		else if( buf[0] == '<')
		{
			if( sscanf(buf+2, "%d %d %*s", &b, &e) != 2 ) {
				fatalf("bad format in %s\n", buf);
			}
		}
		else printf("%s", buf);
	}

	return EXIT_SUCCESS;
}
