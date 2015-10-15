// pull_c -- pull coding sequence from genomic, using CHAP codex file
// starting position: 1 and both ends are included
#include "util.h"
#include "seq.h"

#define AA_LEN_CUTOFF 90
#define AA_CHECK_LEN 300

int b[1000], e[1000], N;

char compl[128];

void find_alternative_start(char *seq, int len, char *flanking_seq, int flanking_len, char dir, int new_b, int new_e, char *chr_name,  char *info);
void find_alternative_stop(char *seq, int len, char dir, int new_b, int new_e, char *chr_name, char *info);
bool is_stop_codon(char *aa);
bool is_start_codon(char *aa);
int check_stop_codon(char *seq, int len, char *next_aa);
int check_start_codon(char *seq);

int main(int argc, char **argv) {
	SEQ *sf;
	uchar *s;
	FILE *fp;
	char buf[1000], chr_name[100], info[1000], dir;
	int i = 0, j = 0, k = 0, n = 0, b = 0, e = 0, B = 0, E = 0;
	int max_len = 0;
	char *cur_seq;
	char *flanking_seq;
	int seq_len = 0, flanking_len = 0;
	int num_nu = 0;
	char next_aa[4];
	int new_b = 0, new_e = 0;
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
	flanking_seq = (char *) ckalloc(sizeof(char) * (max_len+AA_CHECK_LEN+1));
	strcpy(cur_seq, "");
	strcpy(flanking_seq, "");

	while (fgets(buf, 1000, fp)) {
		is_correct_start = false;
		is_correct_stop = false;
		if (sscanf(buf, "%*s %*s %*s %d %d %*s %c %*s", &B, &E, &dir) != 3)
			fatalf("ouch: %s", buf);

		strcpy(cur_seq, "");
		if (dir == '+') {
			n = 0;
			if( B > AA_CHECK_LEN ) b = B - AA_CHECK_LEN;
			else b = 1;

			for ( i = b; i < B; i++ ) {
				flanking_seq[n] = s[i];
				n++;
			}
			flanking_seq[n] = '\0';
			flanking_len = n;

			k = 0;
			for (i = B; i <= E; i++) {
				cur_seq[k] = s[i];
				k++;
			}
			cur_seq[k] = '\0';
			
			if((i + 3) <= seq_len) {
				for( j = 0; j <=2; j++ ) {
					next_aa[j] = s[i];	
					i++;
				}
				next_aa[j] = '\0';
				num_nu = check_stop_codon(cur_seq, k, next_aa);
				if( num_nu != -1 ) {
					is_correct_stop = true;
					new_b = B;
					new_e = E + num_nu;
				}
			} 
		} 
		else {
			if( (E + AA_CHECK_LEN) > seq_len ) e = seq_len;
			else e = E + AA_CHECK_LEN;

			n = 0;	
			for( i = e; i > E; i-- ) {
				flanking_seq[n] = compl[s[i]];
				n++;
			}
			flanking_seq[n] = '\0';
			flanking_len = n;

			k = 0;
			for (i = E; i >= B; i--) {
				cur_seq[k] = compl[s[i]];
				k++;
			}
			cur_seq[k] = '\0';

			if((i - 3) >= 1) {
				for( j = 0; j <=2; j++ ) {
					next_aa[j] = compl[s[i]];	
					i--;	
				}
				next_aa[j] = '\0';
				num_nu = check_stop_codon(cur_seq, k, next_aa);
				if( num_nu != -1 ) {
					is_correct_stop = true;
					new_b = B - num_nu;
					new_e = E;
				}
			} 
		}

		num_nu = check_start_codon(cur_seq);
		if( num_nu != -1 ) {
			is_correct_start = true;
		}

		if( k >= AA_LEN_CUTOFF ) {
			if( is_correct_stop == true ) {
				if( dir == '+' ) {
					find_alternative_start(cur_seq, k, flanking_seq, n, dir, B, new_e, chr_name, info);
				}
				else {
					find_alternative_start(cur_seq, k, flanking_seq, n, dir, new_b, E, chr_name, info);
				}
			}

			if( is_correct_start == true ) {
				if( dir == '+' ) {
					find_alternative_stop(cur_seq, k, dir, new_b, B, chr_name, info);
				}
				else {
					find_alternative_stop(cur_seq, k, dir, E, new_e, chr_name, info);
				}
			}
		}
	}

	free(flanking_seq);
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

int check_stop_codon(char *seq, int len, char *next_aa)
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
	else if( is_stop_codon(next_aa) == true ) {
		res = 3;
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

void find_alternative_start(char *seq, int len, char *flanking_seq, int flanking_len, char dir, int new_b, int new_e, char *chr_name,  char *info)
{
	int i = 0, j = 0;
	bool is_stop = false;
	char prev_aa[4];
	int cur_b = 0, cur_e = 0;
	int loc = 0;
	int e = 0;

	if( (len-AA_LEN_CUTOFF) < 0 ) e = len;
	else e = len-AA_LEN_CUTOFF;

	strcpy(prev_aa, "");

	for( (j = len-1); j >= e; j-- )
	{
		for( i = 2; i >= 0; i-- ) {
			prev_aa[i] = seq[j];	
		}
		if( is_stop_codon(prev_aa) == true ) is_stop = true;
	}

	j = e-1; 
	while(((j-2) >= 0) && (is_stop == false)) {
		for( i = 2; i >= 0; i-- ) {
			prev_aa[i] = seq[j];	
			j--;
		}
		prev_aa[3] = '\0';

		if( is_stop_codon(prev_aa) == true ) is_stop = true;
		else if( is_start_codon(prev_aa) == true ) {
			loc = j+1;
			if( loc != 0 ) {
				if( dir == '+' ) {
					cur_b = new_b + loc;
					cur_e = new_e;
				}
				else {
					cur_b = new_b;
					cur_e = new_e - loc;
				}
				printf("%s maker gene %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
				printf("%s maker CDS %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
			}
		}
	}

	j = flanking_len-1; 
	if( len % 3 == 1 ) {
		if( j >= 1 ) {
			prev_aa[2] = seq[0];
			prev_aa[1] = flanking_seq[j];
			prev_aa[0] = flanking_seq[j-1];
			j = j-2;
		}
	}
	else if( len % 3 == 2 ) {
		if( j >= 0 ) {
			prev_aa[2] = seq[1];
			prev_aa[1] = seq[0];
			prev_aa[0] = flanking_seq[j];
			j = j-1;
		}
	}
	else {
		prev_aa[2] = seq[2];
		prev_aa[1] = seq[1];
		prev_aa[0] = seq[0];
	}
	prev_aa[3] = '\0';

	while(((j-2) >= 0) && (is_stop == false)) {
		if( is_stop_codon(prev_aa) == true ) is_stop = true;
		else if( is_start_codon(prev_aa) == true ) {
//			loc = j+1;
			loc = flanking_len - j;
			if( loc != 0 ) {
				if( dir == '+' ) {
					cur_b = new_b - loc + 1;
					cur_e = new_e;
				}
				else {
					cur_b = new_b;
					cur_e = new_e + loc - 1;
				}
				printf("%s maker gene %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
				printf("%s maker CDS %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
			}
		}

		for( i = 2; i >= 0; i-- ) {
			prev_aa[i] = flanking_seq[j];	
			j--;
		}
		prev_aa[3] = '\0';
	}

}

void find_alternative_stop(char *seq, int len, char dir, int new_b, int new_e, char *chr_name, char *info)
{
	int i = 0, j = len;
	bool is_stop = false;
	char next_aa[4];
	int loc = 0;
	int cur_b = 0, cur_e = 0;

	strcpy(next_aa, "");

	for( j = 0; j < AA_CHECK_LEN; j++ )
	{
		for( i = 0; i <= 2; i++ ) {
			next_aa[i] = seq[j];	
		}
		if( is_stop_codon(next_aa) == true ) is_stop = true;
	}

	j = AA_LEN_CUTOFF; 
	while(((j+2) < len) && (is_stop == false) ) {
		for( i = 0; i < 3; i++ ) {
			next_aa[i] = seq[j];	
			j++;
		}
		next_aa[3] = '\0';
		if( is_stop_codon(next_aa) == true ) is_stop = true;
	}
	
	loc = j-1;

	if( loc < (len-1) ) {
		if( dir == '+' ) {
			cur_b = new_b;
			cur_e = new_b + loc;
		}
		else {
			cur_b = new_e - loc;
			cur_e = new_e;
		}

		if( (loc >= AA_LEN_CUTOFF) && (is_stop == true) ) {
			printf("%s maker gene %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
			printf("%s maker CDS %d %d . %c . %s\n", chr_name, cur_b, cur_e, dir, info);
		}
	}
}
