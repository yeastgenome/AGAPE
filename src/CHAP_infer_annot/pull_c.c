// pull_c -- pull coding sequence from genomic, using CHAP codex file
// starting position: 1 and both ends are included
#include "util.h"
#include "seq.h"

int b[1000], e[1000], N;

char compl[128];

int main(int argc, char **argv) {
	SEQ *sf;
	uchar *s;
	FILE *fp;
	char buf[500], name[100], dir, *status;
	int i = 0, k = 0, B = 0, E = 0, col = 0;

	if (argc != 3)
		fatalf("args: chimp.fa sim4.A=5.out");
	compl['a'] = compl['A'] = 'T';
	compl['c'] = compl['C'] = 'G';
	compl['g'] = compl['G'] = 'C';
	compl['t'] = compl['T'] = 'A';
	sf = seq_get(argv[1]);
	s = SEQ_CHARS(sf) - 1;
	fp = ckopen(argv[2], "r");
	while ((status = fgets(buf, 500, fp)) != NULL && !strchr("<>", buf[0]))
		;
	if (status == NULL)
		fatalf("no genes found in %s", argv[2]);
	while (strchr("<>", dir = buf[0])) {
		if (sscanf(buf+1, "%d %d %s %*s", &B, &E, name) != 3)
			fatalf("ouch: %s", buf);

		printf("> %s %d-%d\n",
	  name, B, E);

		N = 0;
		while ((status = fgets(buf, 500, fp)) != NULL &&
		       !strchr("<>", buf[0]))
			if (sscanf(buf, "%d %d", &(b[N]), &(e[N])) ==  2) {
				if (b[N] > e[N])
					fatalf("exon ends: %d %d\n",
					  b[N], e[N]);
				++N;
			}
		if (dir == '>') {
			for (k = 0; k < N; ++k) {
				if (k > 0) {
					E = e[k-1];
					B = b[k];
//					fprintf(stderr, "intron %c%c .. %c%c\n",
//					  s[E+1], s[E+2], s[B-2], s[B-1]);
				}
				for (col = 0, i = b[k]; i <= e[k]; ++i) {
					putchar(s[i]);
					if (++col == 50) {
						putchar('\n');
						col = 0;
					}
				}
				if (col != 0)
					putchar('\n');
			}
			E = e[N-1];
//			fprintf(stderr, "next codon: %c%c%c (TAA, TAG, TGA?)\n",
//			  s[E+1], s[E+2], s[E+3]);
		} else {
			for (k = N-1; k >= 0; --k) {
				if (k < N-1) {
					E = e[k];
					B = b[k+1];
//					fprintf(stderr, "intron %c%c .. %c%c\n",
//					  compl[s[B-1]], compl[s[B-2]],
//					  compl[s[E+2]], compl[s[E+1]]);
				}
				for (col = 0, i = e[k]; i >= b[k]; --i) {
					putchar(compl[s[i]]);
					if (++col == 50) {
						putchar('\n');
						col = 0;
					}
				}
				if (col != 0)
					putchar('\n');
			}
			B = b[0];
//			fprintf(stderr, "next codon: %c%c%c (TAA, TAG, TGA?)\n",
//			  compl[s[B-1]], compl[s[B-2]], compl[s[B-3]]);
		}
		putchar('\n');
	}

	seq_close(sf);
	fclose(fp);
	return EXIT_SUCCESS;
}
