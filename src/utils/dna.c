/* RCSID("$Id: dna.c,v 1.2 2000/10/31 16:10:07 webb Exp $"); */

/*
* dna -- program to (1) format a DNA sequence, (2) extract a subrange of
* a DNA sequence or (3) form the reverse complement of a DNA sequence.
*
* (1) The command
*	dna filename
* strips out extraneous characters, converts to capital letters, and makes
* lines of 50 nucleotides (except the last line). The input file may begin
* with the character '>', in which case the first line is taken as a comment
* and copied to the output.  Other lines are stripped of their "unknown"
* characters. The "known" characters constitute the string "ABCDGHKMNRSTVWXY".
*
* (2) The command
*	dna from,to filename
* (where "from" and "to" are integers) creates output with the same comment
* line and the indicated subrange of the remaining characters.  The comma can
* be replaced by a string of 0 or more non-digits.  A value 0 for "to" means
* the last position in the sequence (i.e., extract from position "from" onward).
*
* (3) The command
*	dna -c filename
* writes the reverse complement (with the same comment line) on standard out.
*
* The nucleotide ambiguity codes are:
* B = C, G or T
* D = A, G or T
* H = A, C or T
* K = G or T
* M = A or C
* N = A, C, G or T
* R = A or G (purines)
* S = C or G
* V = A, C or G
* W = A or T
* X = A, C, G or T
* Y = C or T (pyrimidines)
*/

#include "util.h"
#include "seq.h"

enum {LINE_LEN=50};

int put_bp(uchar c, int col)
{
	putchar(toupper(c));
	if (++col == LINE_LEN) {
		putchar('\n');
		col = 0;
	}
	return col;
}

void write_seq(uchar *s, int len)
{
	int i, col;

	for (col = i = 0; i < len; ++i, ++s)
		col = put_bp(*s, col);
	if (col != 0)
		putchar('\n');
}

int main(int argc, char *argv[])
{
//	FILE *in;
	SEQ *sf;
	char buf2[1000], *s;
	int from, to;
//	char buf[1000], buf2[1000], *s, *seq, cap, *name, *header;
//	int i, c, len, col, from, to, fake1, fake2;


	if (argc == 2) {					/* Case 1 */
		sf = seq_get(argv[1]);
		if (SEQ_HLEN(sf) != 0)
			printf("%s\n", SEQ_HEAD(sf));
	} else if (argc == 3 && isdigit(argv[1][0])) {		/* Case 2 */
		from = atoi(s=argv[1]);
		while (isdigit(*s))
			++s;
		while (*s && !isdigit(*s))	/* find next digit */
			++s;
		to = atoi(s);
		if (from <= 0 || (to != 0 && to < from))
			fatalf("improper position specification: %s", argv[1]);
		if (to == 0) {
			sf = seq_get(argv[2]);
			to = SEQ_LEN(sf);
			seq_close(sf);
		}
		sprintf(buf2, "%s[%d,%d]", argv[2], from, to);
		sf = seq_get(buf2);
		printf(">%s %d,%d %s\n", argv[0], from, to, argv[2]);
	} else if (argc == 3 && same_string(argv[1], "-c")) {	/* Case 3 */
		sf = seq_get(argv[2]);
		seq_revcomp_inplace(sf);
		printf("> %s %s %s\n", argv[0], argv[1], argv[2]);
	} else {
		fprintf(stderr,
	   "\"%s seq_file\" or \"%s from,to seq_file\" or \"%s -c seq_file\"\n",
	   argv[0], argv[0], argv[0]);
		exit (1);
	}
	write_seq(SEQ_CHARS(sf), SEQ_LEN(sf));
	seq_close(sf);

	return 0;
}
