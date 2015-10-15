#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	int from = 0, to = 0, beg = 0, end = 0;
	int b = 0, e = 0;
	char buf[10000], type[100];
	char cur_buf[10000];
	char scaf_name[100];
	int len = 0, cur_len = 0;
	int strand = 1;
	int num_genes = 0;
	bool is_sname_given = false;
 
	strcpy(scaf_name, "");
	if( argc == 10 ) {
		is_sname_given = true;
		strcpy(scaf_name, argv[9]);
	}
	else if ( argc != 9) {
		fatal("gff2sim4 gff_annotation seq_start seq_end start end gene_name strand seq_len (scaf_name) \n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	
	strcpy(buf, "");
	strcpy(cur_buf, "");
	strcpy(type, "");

	from = atoi(argv[2]);
	to = atoi(argv[3]);
	beg = atoi(argv[4]);
	end = atoi(argv[5]);
	len = atoi(argv[8]);

	cur_len = to - from + 1;
	if(strcmp(argv[7], "+") == 0)  
		strand = 1;
	else 
		strand = 0;
	
	if(beg > end) {
		fatalf("Improper position specification: %d, %d\n", beg, end);
	}

	while(fgets(buf, 10000, f) && (buf[0] != '#')) {};

	while( fgets(buf, 10000, f)) { 
		if((buf[0] == '#') || (buf[0] == '-')) {
		}
		else {
			strcpy(cur_buf, strstr(buf, "AUGUSTUS"));
			if( cur_buf == NULL ) {
				strcpy(cur_buf,  strstr(buf, "augustus"));
				if( cur_buf == NULL ) fatalf("wrong gff format: %s", buf);
			}

			if( sscanf(cur_buf, "%*s %s %d %d %*s", type, &b, &e) != 3 ) {
				fatalf("wrong gff format: %s", cur_buf);
			}
			else if( (strcmp(type, "gene") == 0) ) {

				if( strand == 1 ) {
					if( ((b+from-1) > len) || ((e+from-1) > len) ) {
						fatalf("postions %d or %d exceed the sequence length\n", b+from-1, e+from-1, len);
					}
				}
				else {
					if( ((cur_len-e+from) > len) || ((cur_len-b+from) > len) ) {
						fatalf("postions %d or %d exceed the sequence length\n", b+from-1, e+from-1, len);
					}
				}
				num_genes++;
				if(strand == 1)
					if( is_sname_given == false ) {
						printf("> %d %d %s\n", b + from - 1, e + from - 1, argv[6]);
					}
					else {
						printf("> %d %d %s %s\n", b + from - 1, e + from - 1, argv[6], scaf_name);
					}
			else
					if( is_sname_given == false ) {
						printf("< %d %d %s (complement)\n", cur_len - e + from, cur_len - b + from, argv[6]);
					}
					else {
						printf("< %d %d %s %s (complement)\n", cur_len - e + from, cur_len - b + from, argv[6], scaf_name);
					}
			}
			else if( strcmp(type, "CDS") == 0 ) {
				if( strand == 1 ) {
					printf("%d %d\n", b + from - 1, e + from - 1);
				}
				else {
					printf("%d %d\n", cur_len - e + from, cur_len - b + from);
				}
			}
		}
	}
	
	if( num_genes == 0 ) {
		if( strand == 1 ) {
			if( is_sname_given == false ) {
				printf("> %d %d %s\n", beg, end, argv[6]);
			}
			else {
				printf("> %d %d %s %s\n", beg, end, argv[6], scaf_name);
			}
			printf("%d %d\n", beg, end);
		}
		else {
			if( is_sname_given == false ) {
				printf("< %d %d %s (complement)\n", beg, end, argv[6]);
			}
			else {
				printf("< %d %d %s %s (complement)\n", beg, end, argv[6], scaf_name);
			}
			printf("%d %d\n", beg, end);
		}
	}
	
	fclose(f);

	return EXIT_SUCCESS;
}
