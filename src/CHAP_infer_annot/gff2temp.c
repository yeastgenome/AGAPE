#include "util.h"
#include "util_annot.h"

#define LEN_NAME 100
#define SIMPLE 0
#define SGD 1

int main(int argc, char *argv[])
{
	char buf[5000], cur[5000];
	FILE *f;
	char chr_name[LEN_NAME], dir[LEN_NAME], chr[LEN_NAME], name[LEN_NAME], annot[LEN_NAME], gname[LEN_NAME];
	int b = 0, e = 0;
	int old_b = 0, old_e = 0;
	int num_exons = 0;
	int num_genes = 0;
	int type = SGD;
	
	strcpy(chr, "");
	strcpy(name, "");
	strcpy(annot, "");

	if( argc == 4 ) {
		if( strcmp(argv[3], "SIMPLE") == 0 ) {
			type = SIMPLE;
		} 
		else {
    	fatal("gff2temp gff_file chr_name (SIMPLE)\n");
		}
	}
	else if( argc != 3 ) {
    fatal("gff2temp gff_file chr_name \n");
  }

	strcpy(chr_name, argv[2]);
	f = ckopen(argv[1], "r");

  while(fgets(buf, 5000, f))
  {
    if( (buf[0] == '#') || (buf[0] == '>' ) ) {}
    else if( sscanf(buf, "%s %*s %s %d %d %*s %s %*s %s", chr, annot, &b, &e, dir, cur) != 6 ) {
      fatalf("line in wrong gff format: %s\n", buf);
    }
    else {
			if( strcmp(chr, chr_name) == 0 ) {
				if( strcmp(annot, "gene") == 0 ) {
					if( type == SIMPLE ) strcpy(gname, cur);
					else get_gene_name(cur, gname);
					if( strcmp(dir, "+") == 0 ) {
						if( (num_exons == 0) && (num_genes > 0)) {
							printf("%d %d\n", old_b, old_e);
						}
						printf("> %d %d %s\n", b, e, gname);
						num_exons = 0;
						num_genes++;
						old_b = b;
						old_e = e;
					}
					else if( strcmp(dir, "-") == 0 ) {
						if( num_exons == 0 ) {
							printf("%d %d\n", old_b, old_e);
						}
						printf("< %d %d %s\n", b, e, gname);
						num_exons = 0;
						num_genes++;
						old_b = b;
						old_e = e;
					}
					else {
					}
 	     	}
				else if( strcmp(annot, "CDS") == 0 ) {
					printf("%d %d\n", b, e);
					num_exons++;
				}
 	   	}
		}
  }

	if( (num_genes != 0) && (num_exons == 0) ) {
		printf("%d %d\n", old_b, old_e);
	}
  fclose(f);
	return EXIT_SUCCESS;
}
