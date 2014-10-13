#include "main.h"
#include "util.h"

int main(int argc, char *argv[])
{
	FILE *f, *g;
	char buf[10000];
	char head[MAX_LEN];
	char **names;
	int num_names = 0;
	char cur[MAX_LEN], prev[MAX_LEN];
	int len = 0;
	int i = 0;
	bool *is_in;

	strcpy(buf, "");
	strcpy(head, "");
	strcpy(cur, "");
	strcpy(prev, "");
	if( argc != 3 ) {
		printf("add_unique_only soap_file list_file\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = fopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}

		if(!(g = fopen(argv[2], "r"))) {
			printf("no file %s exists\n", argv[2]);
			return EXIT_FAILURE;
		}
	}

	while(fgets(buf, 10000, f))
	{
		num_names++;
	}

	fseek(f, 0, SEEK_SET);

	i = 0;

	if( num_names > 0 ) {
		names = (char **) ckalloc(num_names * sizeof(char *));
		is_in = (bool *) ckalloc(num_names * sizeof(bool));
	}

	for( i = 0; i < num_names; i++ ) {
		names[i] = (char *) ckalloc(MAX_LEN * sizeof(char));
		strcpy(names[i], "");
		is_in[i] = false;
	}

	i = 0;
	while(fgets(buf, 10000, f))
	{
		if( sscanf(buf, "%s", head) == 1 ) {
			if( i < num_names ) {
				strcpy(names[i], head);
				i++;
			}
			else {
				fatalf("count error: exceeded %d\n", num_names);
			}
		}
	}

	while( fgets(buf, 10000, g) ) {	
		if( sscanf(buf, "%s", head) != 1 ) {
			printf("bad format in %s\n", buf);
			return EXIT_FAILURE;
		}
		else 
		{
			for( i = 0; i < num_names; i++ ) {
				if( strcmp(head, names[i]) == 0 ) {
					is_in[i] = true;
				}
			}
		}
	}

	for( i = 0; i < num_names; i++ ) {
		strcpy(cur, names[i]);
		len = strlen(cur);
		if( names[i][len-2] == '/' ) {
			cur[len-2] = '\0';
		}

		if( (i != 0) && (strcmp(cur, prev) == 0) ) {}
		else {
			if( is_in[i] == false ) {
				if( names[i][len-2] == '/' ) {
					printf("%s/1\n", cur); 	
					printf("%s/2\n", cur); 	
				}
				else {
					printf("%s\n", cur); 	
				}
			}
		}
		strcpy(prev, cur);
	}

	for( i = 0; i < num_names; i++ ) {
		free(names[i]);
	}

	if( num_names > 0 ) {
		free(names);
		free(is_in);
	}

	fclose(g);
	fclose(f);

	return EXIT_SUCCESS;
}
