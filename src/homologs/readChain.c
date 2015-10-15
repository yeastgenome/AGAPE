#include "main.h"
#include "chain.h"
#include "util_i.h"
#include "util_exons.h"
#include "repeats.h"

int main(int argc, char *argv[])
{
	FILE *f;
	struct chain *Chain;
	struct chain *SubChain, *chainToFree;
	struct chain *ch_p, *next_p;
	char buf[NUM_CHARS];
	struct lineFile *lf;
	int i = 0;
	int b = 0, e = 0;
	bool is_null = true;
	struct exons_list *homologs;
	int num_chains = 0;
	int num_homologs = 0;
	struct exons_list *repeats;
	int num_repeats = 0;
	char chr[LEN_NAME];

	strcpy(chr, "");
	if( argc == 3 ) {
		if( (f = ckopen(argv[2], "r")) ) {
			if( fgets(buf, NUM_CHARS, f) ) {
				if( sscanf(buf, "%s %d %d", chr, &b, &e) != 3 ) {
					fatalf("format errors: chr beg end in %s", buf);
				}
			}
			else {
				fatalf("%s is empty\n", argv[2]);
			}
		}
		fclose(f);
	}
	else if( argc != 4 ) {
		fatal("args: chain_file interval_text features_gff_file\n");
	}
	else {
		if( (f = ckopen(argv[2], "r")) ) {
			if( fgets(buf, NUM_CHARS, f) ) {
				if( sscanf(buf, "%s %d %d", chr, &b, &e) != 3 ) {
					fatalf("format errors: chr beg end in %s", buf);
				}
			}
			else {
				fatalf("%s is empty\n", argv[2]);
			}
		}
		fclose(f);
		
		if( (f = ckopen(argv[3], "r")) ) {
			while(fgets(buf, NUM_CHARS, f)) {
				i++;
			}
			num_repeats = i;
			repeats = (struct exons_list *) ckalloc(num_repeats * sizeof(struct exons_list));
			init_exons(repeats, 0, num_repeats-1);	
			fseek(f, 0, SEEK_SET);
			assign_gff_exons_chr(f, repeats, num_repeats, chr);
			quick_sort_inc_exons(repeats, 0, num_repeats-1, POS_BASE);
		}
		else {
			fatalf("file %s invalid\n", argv[4]);
		}
		fclose(f);
	}

	lf = lineFileOpen(argv[1], true);
	Chain = chainRead(lf);
	ch_p = Chain;
	while( (ch_p != NULL) && ((next_p = chainRead(lf)) != NULL) ) {
		ch_p->next = next_p;
		ch_p = ch_p->next;
		i++;
	}

//	printf("Number of chains: %d\n", i);
	i = 0;
	ch_p = Chain;
//	while( (i < NUM_LOOPS) && (ch_p != NULL)  ) {
	while( ch_p != NULL  ) {
//		printf("chain %d: %d-%d\n", ch_p->id, ch_p->tStart, ch_p->tEnd);	
		ch_p = ch_p->next;
		i++;
	}

	num_chains = i;
	homologs = (struct exons_list *) ckalloc(num_chains * sizeof(struct exons_list));
	i = 0;
	f = ckopen(argv[2], "r");
	while( fgets(buf, NUM_CHARS, f) ) { 	
		if( sscanf(buf, "%*s %d %d", &b, &e) != 2 ) {
			fatalf("format errors: chr beg end in %s", buf);
		}
		else {
			ch_p = Chain;

			if( ch_p != NULL ) {
				while( (ch_p != NULL) && (is_null == true) ) {
					chainSubsetOnT(ch_p, b, e, &SubChain, &chainToFree);
					if( SubChain != NULL ) is_null = false;
					ch_p = ch_p->next;
				}
			}

			if( is_null == false ) {
				if( (num_repeats == 0 ) || (is_repeats(repeats, num_repeats, SubChain->tName, SubChain->tStart, SubChain->tEnd) == false) ) {
					homologs[i].reg = assign_I(SubChain->qStart, SubChain->qEnd);
					homologs[i].dir = SubChain->qStrand;
					strcpy(homologs[i].chr, SubChain->qName);
					i++;
				}
//				printf("query: %s %d %d\n", SubChain->qName, SubChain->qStart, SubChain->qEnd);
				if( chainToFree != NULL ) {
					chainFree(&chainToFree);
				}

				while( ch_p != NULL ) {
					chainSubsetOnT(ch_p, b, e, &SubChain, &chainToFree);
					ch_p = ch_p->next;
					if( SubChain != NULL ) {
						if( (num_repeats == 0 ) || ( is_repeats(repeats, num_repeats, SubChain->tName, SubChain->tStart, SubChain->tEnd) == false )) {
							if( SubChain->qStrand == '-' ) {
								homologs[i].reg = assign_I(SubChain->qSize - SubChain->qEnd, SubChain->qSize - SubChain->qStart);
							}
							else {
								homologs[i].reg = assign_I(SubChain->qStart, SubChain->qEnd);
							}
							homologs[i].dir = SubChain->qStrand;
							strcpy(homologs[i].chr, SubChain->qName);
							i++;
						}
//						printf("query: %s %d %d\n", SubChain->qName, SubChain->qStart, SubChain->qEnd);
						if( chainToFree != NULL ) {
							chainFree(&chainToFree);
						}
					}
				}
			}
		}
	}

	num_homologs = i;
	selection_sort_exons(homologs, num_homologs);
//	print_exons_list(homologs, num_homologs);
	num_homologs = remove_redundant_intervals(homologs, num_homologs);
	print_exons_list(homologs, num_homologs);
	free(homologs);
	free(repeats);
	chainFreeList(&Chain);

	fclose(f);
	lineFileClose(&lf);

	return EXIT_SUCCESS;
}
