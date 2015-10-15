#include "main.h"
#include "util.h"

int main(int argc, char *argv[]) {
	FILE *f;
	char buf[10000], type[100];
//	int sum = 0;
	struct n_pair *gnames1, *gnames2;
	int num_gnames1 = 0, num_gnames2 = 0;
	int i = 0, j = 0;
	bool is_in = false;
	bool is_first = false;

	if( argc == 4 ) {
		is_first = true;
	}
	else if( argc != 3 ) {
		fatal("yeast_specific_common two_columns_list single_column_list\n");
	}
	
	strcpy(buf, "");
	strcpy(type, "");

	if((f = ckopen(argv[1], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[1]);
	}
	else {
		while(fgets(buf, 10000, f)) {
			num_gnames1++;
		}
	}

	if( num_gnames1 > 0 ) {
		gnames1 = (struct n_pair *) ckalloc(num_gnames1 * sizeof(struct n_pair));

		for( i = 0; i < num_gnames1; i++ ) {
			strcpy(gnames1[i].name1, "");
			strcpy(gnames1[i].name2, "");
			strcpy(gnames1[i].name2, "");
			gnames1[i].id = 0;
			gnames1[i].len = 0;
		}
	}

	fseek(f, 0, SEEK_SET);
	i = 0;
	while(fgets(buf, 10000, f)) {
		if( sscanf(buf, "%s %s %*s", gnames1[i].name1, gnames1[i].name2) != 2 ) {
			fatalf("wrong format in the gene list: %s", buf);
		}
		i++;
	}
	
	fclose(f);

	if((f = ckopen(argv[2], "r")) == NULL )
	{
		fatalf("Cannot open file %s\n", argv[2]);
	}
	else {
		while(fgets(buf, 10000, f)) {
			num_gnames2++;
		}
	}

	if( num_gnames2 > 0 ) {
		gnames2 = (struct n_pair *) ckalloc(num_gnames2 * sizeof(struct n_pair));

		for( i = 0; i < num_gnames2; i++ ) {
			strcpy(gnames2[i].name1, "");
			strcpy(gnames2[i].name2, "");
			strcpy(gnames2[i].name2, "");
			gnames2[i].id = 0;
			gnames2[i].len = 0;
		}
	}

	fseek(f, 0, SEEK_SET);
	i = 0;
	while(fgets(buf, 10000, f)) {
		if( sscanf(buf, "%s %*s", gnames2[i].name1) != 1 ) {
			fatalf("wrong format in the gene list: %s", buf);
		}
		i++;
	}
	
	fclose(f);

	for( i = 0; i < num_gnames1; i++ ) {
		j = 0; 
		is_in = false;
		while( (j < num_gnames2) && (is_in == false ) ) {
			if( ( (is_first == false ) && (strcmp(gnames1[i].name2, gnames2[j].name1) == 0 )) || ( (is_first == true) && (strcmp(gnames1[i].name1, gnames2[j].name1) == 0 )) ) {
				is_in = true;
			}
			j++;
		}
		if( is_in == false ) {
			printf("%s %s\n", gnames1[i].name1, gnames1[i].name2);
		}
	}

	if( num_gnames1 > 0 ) {
		free(gnames1);
	}

	if( num_gnames2 > 0 ) {
		free(gnames2);
	}

	return EXIT_SUCCESS;
}

