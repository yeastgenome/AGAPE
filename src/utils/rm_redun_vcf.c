#include "main.h"
#include "util.h"

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	char chr[MAX_LEN];
	char cur_chr[MAX_LEN], prev_chr[MAX_LEN];
	int prev_pos = -1, cur_pos = 0;
	bool is_all = false;

	strcpy(buf, "");
	strcpy(chr, "");
	strcpy(prev_chr, "");
	strcpy(cur_chr, "");

	if( argc == 3 ) {
		strcpy(chr, argv[2]);
		if(!(f = ckopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}
	}
	else if( argc != 2 ) {
		printf("rm_redun_vcf vcf_file (or chr)\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = ckopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}
		is_all = true;
	}

	while(fgets(buf, 10000, f))
	{
		if( buf[0] == '#' ) {
			printf("%s", buf);
		}
		else {
			if( sscanf(buf, "%s %d %*s", cur_chr, &cur_pos) != 2 ) {
				printf("wrong args: %s", buf);
				return EXIT_FAILURE;
			}
			else if( ( is_all == true ) || ((is_all == false) && (strcmp(chr, cur_chr) == 0))) 
			{
				if( (strcmp(prev_chr, cur_chr) == 0 ) && (prev_pos == cur_pos ) ) {}
				else {
					printf("%s", buf);
				}
			}
			strcpy(prev_chr, cur_chr);
			prev_pos = cur_pos;
		}
	}

	fclose(f);

	return EXIT_SUCCESS;
}
