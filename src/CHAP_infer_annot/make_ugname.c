#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
	int i = 0, j = 0;
	int res = FALSE;
	char end_three[4];

	len = strlen(name);
	i = len-1;
	
	for( j = 0; j < 3; j++ ) {
		end_three[2-j]=name[i-j];
	}
	end_three[3] = '\0';

	while( (i >= 0) && ( ((name[i] - '0') >= 0) && ((name[i] - '0') <= 9) ) ) i--;

	if( i < 0 ) res = FALSE;
	else if( name[i] == 'P' ) res = TRUE;
	else {
		if( strcmp(end_three, "_ps") == 0 ) {
			res = TRUE;
		}
		else res = FALSE;
	}

	return(res);
}

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[100][500], next[500];
	int b[100], e[100];
	int num_genes = 0;
	int count[100];
	int i = 0, len = 0;
	int is_overlapped = FALSE;
	int is_same = FALSE;
	int valid[100];
	int link[100];
	int num_read = -1;
	int ex_b = 0, ex_e = 1;

	if( argc != 2 ) {
		printf("recal_loc loc_file\n");
		return 1;
	}

	strcpy(buf, "");
	strcpy(next, "");
	for( i = 0; i < 100; i++ ) {
		b[i] = 0;
		e[i] = 1;
		valid[i] = -1;
		link[i] = 0;
		count[i] = 0; 
		strcpy(name[i], "");
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>') 
		{
			fgets(next, 500, f);
			num_read++;
			valid[num_read] = -1;
			if( (sscanf(buf+2, "%d %d %s", &b[num_genes], &e[num_genes], name[num_genes]) != 3 ) || (sscanf(next, "%d %d", &ex_b, &ex_e) != 2)) {
				printf("bad format in %s or %s\n", buf, next);
				return 1;
			}
			else if( (ex_e - ex_b) <= 0 ) {}
			else if( e[num_genes] <= b[num_genes] ) {}
			else if( is_pseudo(name[num_genes]) == TRUE ) {}
			else 
			{
				if( num_genes >= 1 ) {
					is_same = FALSE;
					i = 0;
					while( (i < num_genes) && (is_same == FALSE) )
					{
						if((valid[link[i]] != -1) && (((abs(b[num_genes] - b[i]) <= TH ) && (abs(e[num_genes] - e[i]) <= TH)) || (is_subset(b[num_genes], e[num_genes], b[i], e[i]) == TRUE)) ) is_same = TRUE;
						i++;
					}

					if( is_same == FALSE )
					{
						i = 0;
						while( i < num_genes )
						{
							if( is_subset(b[i], e[i], b[num_genes], e[num_genes]) == TRUE )
							{
								valid[link[i]] = -1;
//					printf("%s: %d-%d is included in %s:%d-%d\n", name[i], b[i], e[i], name[num_genes], b[num_genes], e[num_genes]);
							}
							i++;
						}
					}

					if( is_same == FALSE )
					{
						is_overlapped = FALSE;
						i = 0;
						while( (i < num_genes) && (is_overlapped == FALSE) )  
						{
							if( strcmp(name[num_genes], name[i]) == 0 )
							{
								is_overlapped = TRUE;
								len = strlen(name[num_genes]);
								name[num_genes][len] = 'a' + count[i];
								name[num_genes][len+1] = '\0';
								count[i]++;
							}
							i++;
						}
					}
				}

				if( is_same == FALSE )
				{
//					printf("%c %d %d %s\n", buf[0], b[num_genes], e[num_genes], name[num_genes]);
					valid[num_read] = num_genes;
					link[num_genes] = num_read;
					num_genes++;
					if( num_genes >= 100 ) {
						printf("the number of genes is over 100\n");
						return 1;
					}
				}
			}
		}
		else if(buf[0] == '<')
		{
			fgets(next, 500, f);
			num_read++;
			valid[num_read] = -1;
			if( (sscanf(buf+2, "%d %d %s", &b[num_genes], &e[num_genes], name[num_genes]) != 3 ) || (sscanf(next, "%d %d", &ex_b, &ex_e) != 2)) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else if( e[num_genes] <= b[num_genes] ) {}
			else if( (ex_e - ex_b) <= 0 ) {}
			else if( is_pseudo(name[num_genes]) == TRUE ) {}
			else 
			{
				if( num_genes >= 1 ) {
					is_same = FALSE;
					i = 0;
					while( (i < num_genes) && (is_same == FALSE) )
					{
						if( (valid[link[i]] != -1) && (((abs(b[num_genes] - b[i]) <= TH) && (abs(e[num_genes] - e[i]) <= TH)) || (is_subset(b[num_genes], e[num_genes], b[i], e[i]) == TRUE)) ) is_same = TRUE;
						i++;
					}

					if( is_same == FALSE )
					{
						i = 0;
						while( i < num_genes )
						{
							if( is_subset(b[i], e[i], b[num_genes], e[num_genes]) == TRUE ) {
								valid[link[i]] = -1;
//					printf("%s: %d-%d is included in %s:%d-%d\n", name[i], b[i], e[i], name[num_genes], b[num_genes], e[num_genes]);
							}
							i++;
						}
					}

					if( is_same == FALSE )
					{
						is_overlapped = FALSE;
						i = 0;
						while( (i < num_genes) && (is_overlapped == FALSE) )  
						{
							if( strcmp(name[num_genes], name[i]) == 0 )
							{
								is_overlapped = TRUE;
								len = strlen(name[num_genes]);
								name[num_genes][len] = 'a' + count[i];
								name[num_genes][len+1] = '\0';
								count[i]++;
							}
							i++;
						}
					}
				}

				if( is_same == FALSE )
				{
//					printf("%c %d %d %s (complement)\n", buf[0], b[num_genes], e[num_genes], name[num_genes]);
					valid[num_read] = num_genes;
					link[num_genes] = num_read;
					num_genes++;
					if( num_genes >= 100 ) {
						printf("the number of genes is over 100\n");
						return 1;
					}
				}
			}
		}
//		else if( !is_same ) printf("%s", buf);
	}

	fseek(f, 0, SEEK_SET);
	i = -1;
	while(fgets(buf, 500, f))
	{
		if((buf[0] == '>') || (buf[0] == '<')) {
			i++;	
			if( valid[i] != -1 ) {
				if( buf[0] == '>' ) printf("%c %d %d %s\n", buf[0], b[valid[i]], e[valid[i]], name[valid[i]]);
				else if( buf[0] == '<' ) printf("%c %d %d %s (complement)\n", buf[0], b[valid[i]], e[valid[i]], name[valid[i]]);
			}
		}
		else {
			if( valid[i] != -1 ) printf("%s", buf);
		}
	}

	fclose(f);

	return 0;
}
