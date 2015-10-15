#include "main.h"
#include "util.h"
#include "seq.h"

#define LEN 10000000

char S[LEN];

int main(int argc, char **argv) {
	SEQ *sf;
	int num_len = 0;

	if (argc != 2)
		fatal("args: seq-file");
	
	sf = seq_get(argv[1]);
	strcpy(S, (char *)SEQ_CHARS(sf));
	seq_close(sf);
	num_len = strlen(S);

	printf("%d", num_len);
	return EXIT_SUCCESS;
}
