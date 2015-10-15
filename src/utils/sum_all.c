#include "main.h"
#include "util.h"

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[10000];
	int num_snps = 0, num_pass = 0, num_filter = 0, num_coding = 0, num_coding1 = 0;
	int num_syn = 0, num_non = 0;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0;

	if( argc != 3 ) {
		printf("sum_all snp.summary.txt name\n");
		return EXIT_FAILURE;
	}
	else {
		if(!(f = ckopen(argv[1], "r"))) {
			printf("no file %s exists\n", argv[1]);
			return EXIT_FAILURE;
		}
	}

	while(fgets(buf, 10000, f))
	{
		if( sscanf(buf, "%*s %d %d %d %d %d %*s %*s %d %d", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7 ) {
			fatalf("line in wrong summary file: %s\n", buf);
		}
		else {
			num_snps = num_snps + a1;
			num_pass = num_pass + a2;
			num_filter = num_filter + a3;
			num_coding = num_coding + a4;
			num_coding1 = num_coding1 + a5;
			num_syn = num_syn + a6;
			num_non = num_non + a7;
		}
	}

	printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", argv[2], num_snps, num_pass, num_filter, num_coding, num_coding1, num_syn, num_non);

	fclose(f);

	return EXIT_SUCCESS;
}
