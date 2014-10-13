// pull_c -- pull coding sequence from genomic, using CHAP codex file
// starting position: 1 and both ends are included
#include "util.h"
#include "seq.h"

#define AA_LEN_CUTOFF 90
#define AA_CHECK_LEN 300

int b[1000], e[1000], N;

char compl[128];

bool is_stop_codon(char *aa);
bool is_start_codon(char *aa);
int check_stop_codon(char *seq, int len);
int check_start_codon(char *seq);

int main(int argc, char **argv) {
	SEQ *sf;
	uchar *s;
	FILE *fp;
	char buf[1000], chr_name[100], info[1000], dir;
	int i = 0, k = 0, B = 0, E = 0;
	int max_len = 0;
	char *cur_seq;
	int seq_len = 0;
	int num_nu = 0;
	bool is_correct_stop = false;
	bool is_correct_start = false;

	if (argc != 3)
		fatalf("args: fasta gff");
	compl['a'] = compl['A'] = 'T';
	compl['c'] = compl['C'] = 'G';
	compl['g'] = compl['G'] = 'C';
	compl['t'] = compl['T'] = 'A';
	sf = seq_get(argv[1]);
	s = SEQ_CHARS(sf) - 1;
	seq_len = SEQ_LEN(sf);
	fp = ckopen(argv[2], "r");

	while (fgets(buf, 1000, fp)) {
		if (sscanf(buf, "%s %*s %*s %d %d %*s %c %*s %s", chr_name, &B, &E, &dir, info) != 5) {
			fatalf("ouch: %s", buf);
		}

		if( E > seq_len ) {
			fatalf("gene boundary [%d,%d] over the sequence length %d\n", B, E, seq_len);
		}
		
		if( (E-B+1) > max_len ) max_len = E-B+1;
	}

	fseek(fp, 0, SEEK_SET);
	cur_seq = (char *) ckalloc(sizeof(char) * (max_len+1));
	strcpy(cur_seq, "");

	while (fgets(buf, 1000, fp)) {
		is_correct_start = false;
		is_correct_stop = false;
		if (sscanf(buf, "%*s %*s %*s %d %d %*s %c %*s", &B, &E, &dir) != 3)
			fatalf("ouch: %s", buf);

		strcpy(cur_seq, "");
		if (dir == '+') {
			k = 0;
			for (i = B; i <= E; i++) {
				cur_seq[k] = s[i];
				k++;
			}
			cur_seq[k] = '\0';
			
		} 
		else {
			k = 0;
			for (i = E; i >= B; i--) {
				cur_seq[k] = compl[s[i]];
				k++;
			}
			cur_seq[k] = '\0';
		}

		num_nu = check_stop_codon(cur_seq, k);
		if( num_nu != -1 ) {
			is_correct_stop = true;
		}

		num_nu = check_start_codon(cur_seq);
		if( num_nu != -1 ) {
			is_correct_start = true;
		}

		if( (is_correct_start == true) && (is_correct_stop == true) ) {
			printf("%s", buf);
		}
	}

	free(cur_seq);
	seq_close(sf);
	fclose(fp);
	return EXIT_SUCCESS;
}

int check_start_codon(char *seq)
{
	int res = -1;
	char aa[4];
	strcpy(aa, "");
	aa[0] = seq[0];
	aa[1] = seq[1];
	aa[2] = seq[2];
	aa[3] = '\0';

	if( is_start_codon(aa) == true) {
		res = 0;
	}
	else res = -1;

	return(res);
}

int check_stop_codon(char *seq, int len)
{
	int res = -1;
	char aa[4];
	strcpy(aa, "");
	aa[0] = seq[len-3];
	aa[1] = seq[len-2];
	aa[2] = seq[len-1];
	aa[3] = '\0';

	if( is_stop_codon(aa) == true) {
		res = 0;
	}
	else {
		res = -1;
	}

	return(res);
}

bool is_stop_codon(char *aa)
{
	bool res = false;
  if( (strcmp(aa, "TAA") == 0) || (strcmp(aa, "TAG") == 0) || (strcmp(aa, "TGA") == 0) ) {
		res = true;
	}
	return(res);
}

bool is_start_codon(char *aa)
{
	bool res = false;
	if( strcmp(aa, "ATG") == 0 ) {
		res = true;
	}
	return(res);
}
