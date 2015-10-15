#include "util.h"
#include <stdio.h>
#include <string.h>

int find_gid(char *name, char **gname, int num_genes);
int main(int argc, char *argv[])
{
	FILE *f;
	char buf[1000];
	char **gname;
	char name[500];
	int is_end = 1;
	int num_genes = 0;
	int i = 0;

	strcpy(buf, "");
	if( argc != 3 ) {
		printf("filter_protein_db_fasta fasta gff\n");
		return 1;
	}

	f = fopen(argv[2], "r");
	while(fgets(buf, 1000, f)) num_genes++;

	gname = (char **) ckalloc(sizeof(char *) * num_genes);	
	for( i = 0; i < num_genes; i++ ) {
		gname[i] = (char *) ckalloc(sizeof(char) * 500);
		strcpy(gname[i], "");
	}

	fseek(f, 0, SEEK_SET);
	i = 0;
	while( fgets(buf, 1000, f) ) {
		if( sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %*s %s", gname[i]) != 1 ) {
			printf("wrong gff line: %s\n", buf);
			return 1;
		}
		i++;
	}
	fclose(f);

	f = fopen(argv[1], "r");

	if( fgets(buf, 1000, f) ) is_end = 1;
	else is_end = 0;

	while(is_end != 0)
	{
		if(buf[0] == '>')
		{
			if( sscanf(buf+1, "%s %*s", name) != 1 ) {
				printf("bad format in %s\n", buf);
				return 1;
			}
			else {
				if( find_gid(name, gname, num_genes) != -1 ) {
					printf("%s", buf);
					if(fgets(buf, 1000, f)) is_end = 1;
					else is_end = 0;
					while((is_end != 0 ) &&  (buf[0] != '>')) {
						if( buf[0] != '>' ) printf("%s", buf);

						if( fgets(buf, 1000, f) ) is_end = 1;
						else is_end = 0;
					}
				}
				else {
					if(fgets(buf, 1000, f)) is_end = 1;
					else is_end = 0;
					while((is_end != 0 ) &&  (buf[0] != '>')) {
						if( fgets(buf, 1000, f) ) is_end = 1;
						else is_end = 0;
					}
				}
			}
		}
	}

	fclose(f);
	return 0;
}

int find_gid(char *name, char **gname, int num_genes)
{
	int i = 0;
	int res = -1;

	i = 0;
	while( (i < num_genes) && (strcmp(name, gname[i]) != 0 )) i++;

	if( i >= num_genes ) {
		res = -1;
	}
	else res = i;

 	return(res);
}
