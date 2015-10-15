#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	int from = 0, beg = 0, end = 0;
	int b = 0, e = 0;
	char buf[10000], type[100], seq_len[100];
	char cur_buf[10000];
	int len = 0, cur_len = 0;
	int strand = 1;
	int num_genes = 0;
	bool is_done = false;
 
	if ( argc != 8) {
		fatal("exonerate exonerate_annotation loc_start start end gene_name strand seq_len\n");
	}
	
	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	
	strcpy(buf, "");
	strcpy(cur_buf, "");
	strcpy(type, "");
	strcpy(seq_len, "");

	from = atoi(argv[2]);
	beg = atoi(argv[3]);
	end = atoi(argv[4]);
	len = atoi(argv[7]);

	if(strcmp(argv[6], "+") == 0)  
		strand = 1;
	else 
		strand = 0;
	
	if(beg > end) {
		fatalf("Improper position specification: %d, %d\n", beg, end);
	}

	while(fgets(buf, 10000, f) && (buf[0] != '#')) {};

	while( (is_done == false) && fgets(buf, 10000, f)) { 
		if((buf[0] == '#') || (buf[0] == '-')) {
			if( strstr(buf, "END OF GFF DUMP") != 0 ) {
				is_done = true;
			}
		}
		else {
			if(sscanf(buf, "%s %*s", type) != 1) {
				fatalf("wrong gff format: %s", buf);
			}	
		
			strcpy(cur_buf, strstr(buf, "exonerate"));
			if( cur_buf == NULL ) 
				fatalf("wrong gff format: %s", buf);
		
			if( sscanf(cur_buf, "%*s %s %d %d %s %*s", type, &b, &e, seq_len) != 4 ) {
				fatalf("wrong gff format: %s", buf);
			}
			else if( (strcmp(type, "gene") == 0) ) {

				if( ((b+from) > len) || ((e+from) > len) ) {
					fatalf("postions %d or %d exceed the sequence length %d\n", b+from, e+from, len);
				}

				cur_len = atoi(seq_len);
				num_genes++;
				if(strand == 1)
					printf("> %d %d %s\n", b + from - 1, e + from - 1, argv[5]);
			else
				printf("< %d %d %s (complement)\n", cur_len - e + 1 + from, cur_len - b + 1 + from, argv[5]);
			}
			else if( strcmp(type, "exon") == 0 ) {
				if( strand == 1 ) {
					printf("%d %d\n", b + from - 1, e + from - 1);
				}
				else {
					printf("%d %d\n", cur_len - e + 1 + from, cur_len - b + 1 + from);
				}
			}
		}
	}
	
	if( num_genes == 0 ) {
		if( strand == 1 ) {
			printf("> %d %d %s\n", beg, end, argv[5]);
			printf("%d %d\n", beg, end);
		}
		else {
			printf("< %d %d %s (complement)\n", beg, end, argv[5]);
			printf("%d %d\n", beg, end);
		}
	}
	
	fclose(f);

	return EXIT_SUCCESS;
}
