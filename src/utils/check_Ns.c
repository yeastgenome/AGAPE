#include "main.h"
#include "util.h"
#include "seq.h"

#define LEN 10000000
#define LONG_LEN 300

char S[LEN];

int main(int argc, char **argv) {
	SEQ *sf;
	int num_len = 0;
//	int col = 0;
	int i = 0;
	int N;
	int num_N = 0;
	int num_no_N = 0;
	bool is_N = false;
	int num_valid_chunks = 0;
	bool is_valid = true;
	float ratio = (float) 0;
//	int len = 0;

	if (argc != 3)
		fatal("args: seq-file ratio");
	
	ratio = atof(argv[2]);
	sf = seq_get(argv[1]);
	N = SEQ_LEN(sf);
	strcpy(S, (char *)SEQ_CHARS(sf));
	seq_close(sf);
	num_len = strlen(S);
	for( i = 0; i < num_len; i++ ) {
		if( strchr("NBDHKMRSVWY", toupper(S[i])) ) {
			num_N++;
			if( is_N == false ) {
				if( num_no_N >= LONG_LEN ) {
					num_valid_chunks++;
				}
				num_no_N = 0;
				is_N = true;
			}
		}
		else {
			is_N = false;
			num_no_N++;
		}
	}

	if( is_N == false ) {
		if( num_no_N >= LONG_LEN ) num_valid_chunks++;
	}

	if( num_valid_chunks == 0 ) is_valid = false;
	if( ((float)num_N / (float)num_len) >= ratio ) is_valid = false;

	if( is_valid == false ) printf("invalid"); 
	else printf("valid");
	return EXIT_SUCCESS;
}
