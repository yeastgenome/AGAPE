#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"

#define BIG 1000000
#define NUM_COL 13

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[BIG];
	int **num;
	int num_rows = 0, i = 0, j = 0;
	int *sum;

	if( argc != 2 ) {
		printf("args: summary_file\n");
		return EXIT_FAILURE;
	}
	f = fopen(argv[1], "r");
		
	while(fgets(buf, BIG, f))
	{
		i++;
	}
	num_rows = i;

	num = (int **) ckalloc(NUM_COL * sizeof(int *));
	for( i = 0; i < NUM_COL; i++ ) {
		num[i] = (int *) ckalloc(num_rows * sizeof(int));
		for( j = 0; j < num_rows; j++ ) {
			num[i][j] = 0;
		}
	}

	sum = (int *) ckalloc(NUM_COL * sizeof(int));
	for( j = 0; j < NUM_COL; j++ ) sum[j] = 0;

	fseek(f, 0, SEEK_SET);

	i = 0;
	while(fgets(buf, BIG, f))
	{
		if( sscanf(buf, "%*s %d %d %d %d %d %d %d %d %d %d %d %d %d", &num[0][i], &num[1][i], &num[2][i], &num[3][i], &num[4][i], &num[5][i], &num[6][i], &num[7][i], &num[8][i], &num[9][i], &num[10][i], &num[11][i], &num[12][i]) != 13 ) {
			printf("bad format in %s\n", buf);
			return EXIT_FAILURE;
		}
		else {
			i++;
		}
	}

	if( i > num_rows ) {
		printf("counting errors: %d vs %d\n", i, num_rows);
		return EXIT_FAILURE;
	}

	for( i = 0; i < NUM_COL; i++ ) {
		for( j = 0; j < num_rows; j++ ) {
			sum[i] = sum[i] + num[i][j];
		}
	}

	printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", sum[0], sum[1], sum[2],sum[3], sum[4], sum[5], sum[6], sum[7], sum[8], sum[9], sum[10], sum[11], sum[12]);           

	for( i = 0; i < NUM_COL; i++ ) free(num[i]);
	free(num);
	free(sum);

	fclose(f);

	return EXIT_SUCCESS;
}
