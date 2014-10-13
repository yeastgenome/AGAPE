#include "main.h"
#include "util.h"

bool is_in_glist(char *gname, struct n_pair *glist, int num_genes);

int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000];
	char column[1000];
	char column1[1000], column2[1000];
	struct n_pair *gnames;
	int num_gnames = 0;
	char name1[100], species1[100], name2[100], species2[100];
	float pid1 = (float) 0, pid2 = (float) 0; 
	int i = 0, j = 0, cur_id = 0;
	bool is_in = false;
	bool is_first = false;
	char sp_name[100];

	if( argc != 3 ) {
		fatal("get_glsit gff sp_name\n");
	}

	strcpy(sp_name, argv[2]);
	
	strcpy(buf, "");
	strcpy(column, "");
	strcpy(column1, "");
	strcpy(column2, "");
	strcpy(name1, "");
	strcpy(name2, "");
	strcpy(species1, "");
	strcpy(species2, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		if( strstr(buf, ";") != NULL ) num_gnames++;
		num_gnames++;
	}

	if( num_gnames > 0 ) {
		gnames = (struct n_pair *) ckalloc(num_gnames * sizeof(struct n_pair));

		for( i = 0; i < num_gnames; i++ ) {
			strcpy(gnames[i].name1, "");
			strcpy(gnames[i].name2, "");
			gnames[i].id = 0;
			gnames[i].len = 0;
		}
	}

	i = 0;
	fseek(f, 0, SEEK_SET);
	while(fgets(buf, 10000, f)) {
		if( sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %*s %s", column) != 1 ) {
			fatalf("wrong format in gff: %s", buf);
		}
		else {
			if( strstr(buf, ";") != NULL ) {
				if( sscanf(column, "%[^;];%s", column1, column2) != 2 ) {
					fatalf("wrong format in gff: %s", column);
				}
				else {
					if( sscanf(column1, "%[^,],%[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%f", name1, species1, pid1) != 3 ) {
						fatalf("wrong format in gff: %s", column);
					}
				
					if( sscanf(column2, "%[^,],%[^,],%*[^,],%*[^,],%*[^,],%*[^,],%*[^,],%f", name2, species2, pid2) != 3 ) {
						fatalf("wrong format in gff: %s", column);
					}

					if( (pid2 - pid1) <= (float)2 ) {
						if( (strcmp(species2, sp_name) == 0 ) && (is_in_glist(name2, species2, gnames, num_gnames) == false )) {
							strcpy(gnames[i].name1, name2);
							strcpy(gnames[i].name2, species2);
							gnames[i].len = (int)pid2;
							i++;
						}
					}

					if( (strcmp(species1, sp_name) == 0 ) && (is_in_glist(name1, species1, gnames, num_gnames) == false )) {
						strcpy(gnames[i].name1, name1);
						strcpy(gnames[i].name2, species1);
						gnames[i].len = (int)pid1;
						i++;
					}
				}
			}
			else {
				if( (strcmp(species1, sp_name) == 0 ) && (is_in_glist(name1, gnames, num_gnames) == false )) {
					strcpy(gnames[i].name1, name1);
					strcpy(gnames[i].name2, species1);
					gnames[i].len = (int)pid1;
					i++;
				}
			}
		}	
	}

	fclose(f);

	for( i = 0; i < num_gnames1; i++ ) {
		printf("%s\n", gnames[i].name);
	}

	if( num_gnames > 0 ) {
		free(gnames);
	}

	return EXIT_SUCCESS;
}

bool is_in_glist(char *gname, struct n_pair *glist, int num_genes)
{
	int i = 0;
	bool res = false;

	while( ( i < num_genes) && (res == false) ) {
		if( strcmp(gname, glist[i].name1) == 0 ) {
			res = true;
		}
		i++;
	}
}
