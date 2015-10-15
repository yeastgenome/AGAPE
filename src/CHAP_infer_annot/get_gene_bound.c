#include <stdio.h>
#include <string.h>
#include <math.h>

#define TRUE 1
#define FALSE 0
#define TH 30

int is_subset(int b1, int e1, int b2, int e2)
{
	if( ((b2 - b1) < TH)  && ((e1 - e2) < TH) ) return(TRUE);
	else return(FALSE);
}

int is_pseudo(char *name)
{
	int len = 0;
	int i = 0;
	int res = FALSE;

	len = strlen(name);
	i = len-1;
	while( (i >= 0) && ( ((name[i] - '0') >= 0) && ((name[i] - '0') <= 9) ) ) i--;

	if( i < 0 ) res = FALSE;
	else if( name[i] == 'P' ) res = TRUE;
	else res = FALSE;

	return(res);
}

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500];
	int b = 0, e = 1;

	strcpy(buf, "");
	strcpy(name, "");

	if( argc != 2 ) {
		printf("get_gene_bound loc_file\n");
		return 1;
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if((buf[0] == '>') || (buf[0] == '<')) {
			sscanf(buf+2, "%d %d %s %*s", &b, &e, name);
			if( buf[0] == '>' ) printf("%c %d %d %s\n", buf[0], b, e, name);
			else if( buf[0] == '<' ) printf("%c %d %d %s\n", buf[0], b, e, name);
			printf("%d %d\n", b, e);
		}
	}

	fclose(f);

	return 0;
}
