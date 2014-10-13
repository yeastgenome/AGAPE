#include "util.h"

void print_exons_reverse(int **exons, int num_exons);
int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int **exons;
	int num_exons = 0, max_num_exons = 0;
	int i = 0;
	bool is_prev_reverse = false;

	strcpy(buf, "");
	if( argc != 2 ) {
		printf("reverse_exon_order loc_file\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = fopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}
	}

	max_num_exons = 0;
	while(fgets(buf, 10000, f))
	{
		if( (buf[0] == '>') || (buf[0] == '<') ) {
			if( max_num_exons < num_exons ) max_num_exons = num_exons;
			num_exons = 0;
		}
		else {
			num_exons++;
		}	
	}
	if( max_num_exons < num_exons ) max_num_exons = num_exons;

	fseek(f, 0, SEEK_SET);

	if( max_num_exons > 0 ) {
		exons = (int **) ckalloc(max_num_exons * sizeof(int *));
		for( i = 0; i < max_num_exons; i++ ) {
			exons[i] = (int *) ckalloc(2 * sizeof(int));
		}
	}

	i = 0;

	while(fgets(buf, 10000, f))
	{
		if(buf[0] == '<') {
			if( is_prev_reverse == true ) print_exons_reverse(exons, num_exons);
			printf("%s", buf);
			num_exons = 0;
			is_prev_reverse = true;
		}
		else if( buf[0] == '>' ) {
			if( is_prev_reverse == true ) print_exons_reverse(exons, num_exons);
			is_prev_reverse = false;
			printf("%s", buf);
			num_exons = 0;
		}
		else {
			if( is_prev_reverse == true ) {
				if( num_exons >= max_num_exons ) {
					fatalf("error: %d exceeds max_num\n", num_exons);
				}

				if( sscanf(buf, "%d %d %*s", &exons[num_exons][0], &exons[num_exons][1]) != 2 ) {
					fatalf("format error: %s", buf);
				}
				num_exons++;
			}
			else {
				printf("%s", buf);
			}
		}	
	}
	if( is_prev_reverse == true ) print_exons_reverse(exons, num_exons);

	for( i = 0; i < max_num_exons; i++ ) {
		free(exons[i]);
	}

	if( max_num_exons > 0 ) {
		free(exons);
	}

	fclose(f);

	return EXIT_SUCCESS;
}

void print_exons_reverse(int **exons, int num_exons)
{
	int i = 0;

	if( exons[num_exons-1][1] < exons[0][0] ) {
		for( i = num_exons-1; i >= 0; i-- ) {
			printf("%d %d\n", exons[i][0], exons[i][1]);
		}
	}
	else {
		for( i = 0; i < num_exons; i++ ) {
			printf("%d %d\n", exons[i][0], exons[i][1]);
		}
	}
}
