#include "util.h"
#include "util_i.h"

#define LEN_CUTOFF_RATIO 0.7

int main(int argc, char *argv[]) {
	FILE *f;
	int b = 0, e = 0;
	struct I cur_reg, query, common_reg;
	char query_name[500];
	char buf[10000], type[100];
	char scaf_name[500];
	float len_cutoff_ratio = LEN_CUTOFF_RATIO;
	char *is_null = NULL;
	char gname[500];
	char sign = '+';

	strcpy(query_name, "");
	strcpy(scaf_name, "");
	strcpy(gname, "");
	strcpy(type, "");
	strcpy(buf, "");
 
	cur_reg = assign_I(0, 1);
	query = assign_I(0, 1);
	common_reg = assign_I(0, 1);

	if( argc == 7 ) {
		len_cutoff_ratio = atof(argv[5]);
		strcpy(gname, argv[6]);
	}
	else if( argc == 6 ) {
		len_cutoff_ratio = atof(argv[5]);
	}
	else if ( argc != 5) {
		fatal("gff_subset gff_annotation scaf_name b e (ratio)\n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}

	b = atoi(argv[3]);
	e = atoi(argv[4]);
	strcpy(query_name, argv[2]);
	query = assign_I(b, e);
	is_null = fgets(buf, 10000, f);
	while( is_null != NULL ) { 
		if( sscanf(buf, "%s %*s %s %d %d %*s %c %*s", scaf_name, type, &b, &e, &sign) != 5 ) {
			fatalf("wrong gff line: %s", buf);
		}
		else
		{
			cur_reg = assign_I(b, e);
			if( (strcmp(type, "gene") == 0 ) && (strcmp(scaf_name, query_name) == 0) ) 
			{
				if( proper_overlap(cur_reg, query) == true ) {
					common_reg = intersect(cur_reg, query);
					if( width(common_reg) >= (int)(len_cutoff_ratio * ((float)(width(cur_reg)) )) ) {
						if( argc == 7 ) {
							printf("%s\tmaker\t%s\t%d\t%d\t.\t%c\t.\t%s\n", scaf_name, type, b, e, sign, gname);
						}
						else {
							printf("%s", buf);
						}
						is_null = fgets(buf, 10000, f);
						if( sscanf(buf, "%s %*s %s %d %d %*s", scaf_name, type, &b, &e) != 4 ) 
						{
							fatalf("wrong gff line: %s", buf);
						}
						else { 
							while( (is_null != NULL ) && (strcmp(type, "CDS") == 0) ) {
								if( argc == 7 ) {
									printf("%s\tmaker\t%s\t%d\t%d\t.\t%c\t.\t%s\n", scaf_name, type, b, e, sign, gname);
								}
								else {
									printf("%s", buf);
								}
								is_null = fgets(buf, 10000, f);
								if( sscanf(buf, "%s %*s %s %d %d %*s", scaf_name, type, &b, &e) != 4 ) 
								{
									fatalf("wrong gff line: %s", buf);
								}
							}
						}
					}
					else {
						is_null = fgets(buf, 10000, f);
					}
				}
				else {
					is_null = fgets(buf, 10000, f);
				}
			}	
			else {
				is_null = fgets(buf, 10000, f);
			}
		}
	}
	
	fclose(f);

	return EXIT_SUCCESS;
}
