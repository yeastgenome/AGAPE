#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "util.h"

#define TRUE 1 
#define FALSE 0

#define NORMAL 0
#define DEBUG_MODE 1

void get_one_aa(int b, char *fname, char *res);

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500], name[500];
	int b = 0, e = 1;
	char beg = 0, end = 1;
	int len = 0;
	char ext_b[4], ext_e[4];
	int total_len = 0;
	int org_len = 0;
	int debug_mode = NORMAL;

	strcpy(buf, "");
	strcpy(name, "");

	if( argc == 5 ) debug_mode = DEBUG_MODE;
	else if( argc != 4 ) {
		fatal("filter_out seq_file len\n");
	}

	org_len = atoi(argv[3]);

	f = fopen(argv[1], "r");

	fgets(buf, 500, f);
	if(buf[0] == '>') 
	{
		if( sscanf(buf+2, "%s %d-%d", name, &b, &e) != 3 )
		{
			fatalf("bad format in %s\n", buf); 
		}
		else
		{
			fgets(buf, 500, f);
			beg = buf[0];
		}
	}

	while(fgets(buf, 500, f)) total_len = total_len + strlen(buf) - 1;

	len = strlen(buf);
	end = buf[len-2];

	if( (org_len >= 300) && (total_len < (org_len/3) ) ) printf("f");
	else {
		if( (beg == 'M') && (end == '*') ) { printf("t"); }
		else if( beg == 'M') {
			get_one_aa(e+1, argv[2], ext_e);
			if( debug_mode == DEBUG_MODE ) printf("%s\n", ext_e);
			if( (strcmp(ext_e, "TAG") == 0) || (strcmp(ext_e, "TGA") == 0) || (strcmp(ext_e, "TAA") == 0) ) printf("P");
			else printf("f");
		}
		else if( end == '*' ) {
			get_one_aa(b-3, argv[2], ext_b);
			if( debug_mode == DEBUG_MODE ) printf("%s\n", ext_e);
			if( strcmp(ext_b, "ATG") == 0 ) printf("M");
			else printf("f");
		}
		else {
			get_one_aa(b-3, argv[2], ext_b);
			get_one_aa(e+1, argv[2], ext_e);
			if( ((strcmp(ext_e, "TAG") == 0) || (strcmp(ext_e, "TGA") == 0) || (strcmp(ext_e, "TAA") == 0)) && (strcmp(ext_b, "ATG") == 0 ) ) printf("b");
			else printf("f");
		}
	}

	return EXIT_SUCCESS;
}

void get_one_aa(int b, char *fname, char *res)
{
	FILE *seq;
	char buf[500];
	int old_count = 0, count = 0;
	int len = 0, l = 0, i = 0, j = 0;

	seq = fopen(fname, "r");

	fgets(buf, 500, seq);
	while((count < b) && fgets(buf, 500, seq)) 
	{
		old_count = count;
		len = strlen(buf);
		count = count + len - 1;		
	}

	l = b - old_count - 1;

	if( (len - 1 - l) >= 3 ) {
		for( i = 0; i < 3; i++ ) res[i] = toupper(buf[l + i]);
	}
	else if( ((len - 1 - l) < 3) && (( len - 1 - l ) > 0 )) {
		for( i = 0; i < (len - 1 - l); i++ ) res[i] = toupper(buf[l + i]);

		fgets(buf, 500, seq);

		for( j = 0; j < (3 - (len - 1 - l)); j++ ) res[i+j] = toupper(buf[j]);
	}
	else {
		strcpy(res, "");
	}
	res[3] = '\0';
}
