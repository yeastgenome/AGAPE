#include "main.h"
#include "util.h"

void print_genes_in_scaf(FILE *f, char *seq_name);

int main(int argc, char *argv[]) {
	FILE *f;

	if( argc != 3) {
		fatal("extract_scaf_codex codex seq_name\n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	print_genes_in_scaf(f, argv[2]);

	fclose(f);
	return EXIT_SUCCESS;
}

void print_genes_in_scaf(FILE *f, char *seq_name)
{
	char buf[10000], name[100], scf_name[100];
	int num_genes = 0;
	int cur_exons_count = 0;
	int b = 0, e = 0;
	bool is_in = false;

	while(fgets(buf, 10000, f))
	{
		if( buf[0] == '#' ) {}
		else if( (buf[0] == '>') || (buf[0] == '<') )
		{
			if( ( num_genes > 0 ) && ( cur_exons_count == 0 ) ) {
				printf("%d %d\n", b, e);	
			}

			if( sscanf(buf, "%*s %d %d %s %s %*s", &b, &e, name, scf_name) == 4 ) {
				if( strcmp(scf_name, seq_name) == 0 ) {
					printf("%s", buf);
					is_in = true;
				}
				else is_in = false;
				cur_exons_count = 0;
				num_genes++;
			}
			else {
				printf("seq name not given in %s\n", buf);	
			}
		}
		else {
			if( is_in == true ) {
				printf("%s", buf);
			}
			cur_exons_count++;
		}
	}

	if( ( num_genes > 0 ) && ( cur_exons_count == 0 ) ) {
		printf("%d %d\n", b, e);	
	}
}
