#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	int b = 0, e = 0;
	char buf[10000], type[100];
	char column[10000];
	char rest[10000];
	char scaf_name[100];
	char strand = '+';
	bool is_cds = false;
	bool is_orf = false;
	char item1[500], item2[500];
	char gname[2000];
	char source[100];
	char name1[500], name2[500];
	bool is_name_given = false;

	strcpy(item1, "");
	strcpy(item2, "");
	strcpy(column, "");
	strcpy(gname, "");
	strcpy(source, "");
	strcpy(name1, "");
	strcpy(name2, "");
 
	if( argc == 4 ) {
		is_name_given = true;
		if( strcmp(argv[2], "CDS") == 0 ) {
			is_cds = true;
		}
		else if( strcmp(argv[2], "ORF") == 0 ) {
			is_orf = true;
		}
		else {
			fatal("gff2codex gff_annotation (CDS|ORF)\n");
		}
	}
	else if ( argc == 3 ) {
		if( strcmp(argv[2], "CDS") == 0 ) {
			is_cds = true;
		}
		else if( strcmp(argv[2], "ORF") == 0 ) {
			is_orf = true;
		}
		else {
			fatal("gff2codex gff_annotation (CDS|ORF)\n");
		}
	}
	else if ( argc != 2) {
		fatal("gff2codex gff_annotation (CDS)\n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	
	strcpy(buf, "");

	while( fgets(buf, 10000, f)) { 
		strcpy(gname, "");
		strcpy(column, "");
		if((buf[0] == '#') || (buf[0] == '-')) {
		}
		else {
			if( sscanf(buf, "%s %s %s %d %d %*s %c %*s %s %s %*s", scaf_name, source, type, &b, &e, &strand, column, rest) == 8 ) {
			}
			
			if( sscanf(buf, "%s %s %s %d %d %*s %c %*s %s %*s", scaf_name, source, type, &b, &e, &strand, column) != 7 ) {
				fatalf("wrong gff format: %s", buf);
			}
			else if( strcmp(type, "gene") == 0) 
			{
				strcpy(name1, "UNDEF");
				if(strcmp(source, "SGD") == 0) {
					if( (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%s", item1, item2) != 2) ) {
						fatalf("wrong SGD gff column: %s", column);
					}
					else 
					{
						if( sscanf(item1, "%*[^=]=%s", name1) != 1 ) {
              fatalf("wrong gff column: %s", column);
            }
						sprintf(gname, "%s;%s", item1, item2);
					}
				}
				else strcpy(gname, column);

				if(strand == '+') {
					if( is_name_given == false ) {
						printf("> %d %d %s %s\n", b, e, gname, scaf_name);
					}
					else {
						printf("> %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					if( is_cds == false ) printf("%d %d\n", b, e);
				}
				else {
					if( is_name_given == false) {
						printf("< %d %d %s %s\n", b, e, gname, scaf_name);
					}
					else {
						printf("< %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					if( is_cds == false ) printf("%d %d\n", b, e);
				}
			}
			else if((strcmp(type, "match") == 0) || (strcmp(type, "match_part") == 0) ) {
				strcpy(name1, "UNDEF");
				if(strcmp(source, "SGD") == 0) {
					if( (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%s", item1, item2) != 2) ) {
						fatalf("wrong SGD gff column: %s", column);
					}
					else {
						if( sscanf(item1, "%*[^=]=%s", name1) != 1 ) {
              fatalf("wrong gff column: %s", column);
            }
						sprintf(gname, "%s;%s", item1, item2);
					}
				}
				else strcpy(gname, column);

				if(strand == '+') {
					if( is_name_given == false ) {
						printf("> %d %d %s %s\n", b, e, gname, scaf_name);
					}
					else {
						printf("> %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					printf("%d %d\n", b, e);
				}
				else {
					if( is_name_given == false ) {
						printf("< %d %d %s %s\n", b, e, gname, scaf_name);
					}
					else {
						printf("< %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					printf("%d %d\n", b, e);
				}
			}
			else if( strcmp(type, "CDS") == 0) {
				strcpy(name2, "UNDEF");
				if(strcmp(source, "SGD") == 0) {
					if( (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%*[^;];%s", item1, item2) != 2) && (sscanf(column, "%[^;];%*[^;];%s", item1, item2) != 2) ) {
						fatalf("wrong SGD gff column: %s", column);
					}
					else { 
						if( sscanf(item1, "%*[^=]=%s", item2) != 1 ) {
              fatalf("wrong gff column: %s", column);
            }
            else {
					  	if( sscanf(item2, "%[^_]_%*s", name2) != 1 )
              	fatalf("wrong gff column: %s", column);
						}
					}
				}

				if( (is_cds == true) && (strcmp(name1, name2) == 0) ) printf("%d %d\n", b, e);
			}
			else if( (is_orf == true ) && (strcmp(type, "ORF") == 0) ) {
				if(strand == '+') {
					if( is_name_given == false ) {
						printf("> %d %d %s %s\n", b, e, rest, scaf_name);
					}
					else {
						printf("> %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					printf("%d %d\n", b, e);
				}
				else {
					if( is_name_given == false ) {
						printf("< %d %d %s %s\n", b, e, rest, scaf_name);
					}
					else {
						printf("< %d %d %s %s\n", b, e, argv[3], scaf_name);
					}
					printf("%d %d\n", b, e);
				}
			}
		}
	}
	
	fclose(f);

	return EXIT_SUCCESS;
}
