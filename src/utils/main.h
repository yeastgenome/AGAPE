#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include "util.h"

/* If local alignments are expected to be very strong, i.e. most of
 * alignments have high similarities, STRICT is set as a low number like 2
 * otherwise, STRICT is set as a high number 
 * ALT_EFFEC_VALUE is also a threshold value that should be changed according
 * to the type of local alignments.
 * If it has alignments in low percentage identities, it is 300. // for loose alignments
 * Otherwise, ALT_EFFEC_VALUE is set to be zero.
 */

#define SGD 1
#define ENSEMBL 2
#define MAKER 3
#define MULTI_CDS 4

// #define STRICT 2 //for strong local alignments
#define STRICT 8 // for loose local alignments
// #define LOOSE 8 // for strong local alignments
#define LOOSE 10 
#define LOOSE_RUN 5

// running mode
#define INF_DUP 0 // generating only duplication events
#define ONE_TO_ONE 1 // generating one-to-one mappings
#define MANY_TO_MANY 2 // generating many-to-many mappings
#define ONE_TO_MANY 3 // generating one-to-many mappings
#define INC_CONV 4 // including conversions

#define ALT_EFFEC_VALUE 0 // for strong local alignments
// #define ALT_EFFEC_VALUE 300 // for loose local alignments

#define MIN_GAP 100
#define C_OP_TH 4000 // Threshold value of the distance between two close overlapped alignments
#define WIN_SIZE 100 // a window size of comparing the similarity level between two alignments
#define BIG 1000000
#define PID_TH 50 // a threshold of a percentage identity in a self-alignment 

#define SM_TH 30
#define FIRST_SEQ 1
#define SECOND_SEQ 2
#define FWD 0
#define BWD 1

#define CANDI_MODE 0
#define ORTHO_MODE 1
#define BACK_MODE 2

#define MAX_DIFF 100
#define MIN_DIFF -100
#define PID_DIFF_TH 5

#define DEALLOC 1
#define NO_DEALLOC 2

#define BEFORE_SP 1
#define AFTER_SP 2

#define TRUE 1
#define FALSE 0

#define READY 2
 
#define MCMC 0
#define GREEDY 1

#define NUM_SP 3 // the number of species including an out-group species
#define MAX_NUM_GENES 1000 // the maximum number of genes in one tandem gene cluster
#define MAX_NEWICK 10000 // the max length of characters in a newick format tree
#define MAX_LEN 2000 // the max length of a line in a gene annotation file
#define LEN_NAME 100 // the length of a gene name
#define LEN_LONG_NAME 10000 // the length of a gene name
#define STACK_SIZE 100
#define MAX_NAME 100
#define MAX_OPS 500 // the max number of operations for one tandem gene cluster
#define NUM_GENES 200 // the max number of genes which are involved in a gene cluster
#define INC_TH 0.9 // the threshold of percentage to decide the inclusion of a region for other region
#define TH_DIS 500
#define TH_LEN 200

#define F_PAR 0
#define F_LEAF 1
#define VISITED -1

#define DP 0 // duplication
#define CV 1 // conversion
#define SP 2 // speciation
#define SP_W_LOSS 3 // speciation with loss
#define DP_W_LOSS 4 // duplication with loss
#define DP_OF_LEAVES 5 // duplication between closest two species
#define LEAF 6 // leaf node
#define ROOT 7 // root node
#define DONE 8 // node done temporarily
#define DP_CON 9 // confirmed duplication

//enum bool_type {false, true};
//typedef enum bool_type bool;

#define BASE_NUM 20

#define SELF1 0
#define SELF2 1
#define PAIR 2
#define SELF 3
#define INIT_SELF1 10
#define INIT_SELF2 11

#define DEL_EXIST -2
#define NO_EXIST -1
#define SP_EVENT -3

#define X_VAL 0
#define Y_VAL 1

#define MAX_L_ALG 1000
#define MAX_NUM_ID 1000
#define LG_LEN 3000000
#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

#define W_SID 0
#define W_FID 1
#define H_SID 2
#define H_FID 3
#define OW_SID 4
#define OH_FID 5

#define MBP 150

#define FIRST_RUN 1
#define SECOND_RUN 2

#define RM_TH_UNIT 5

#define SP_1 1
#define SP_2 2

#define DELETED 2 // deleted alignment
#define TEMP_HIDDEN 30 // temporally hidden 
#define TEMP_HIDDEN_COMP 31 // temporally hidden complement
#define NON_ORTHO 10
#define OUT_PAR 11
#define OUT_PAR_REV 12
#define ORTHO 20
#define ORTHO_COMP 21
#define ORTHO_CANDI 23
#define ORTHO_COMP_CANDI 24
#define IN_PAR 13
#define IN_PAR_REV 14
#define IN_PAR_IN_X 15
#define IN_PAR_IN_Y 16

#define COPY_OWN 1
#define COPY_OWN_INV 2
#define NON_COPY -1

#define LEFT_LOCK 1
#define RIGHT_LOCK 2
#define BOTH_LOCK 3

#define PERFECT_ID 99 // the basis of the perfect match ID
#define FRAC_P_MATCHES 0.90 // the fraction of the perfect matches

#define MIN_LEN_FOR_RECAL_PID 200
#define Num_Ids 10
// #define NUM_RUN 2 // the number of trial to change to loosen the threshold value
// #define NUM_RUN 3 // the number of trial to change to loosen the threshold value
#define Len_Seq 200000 // the length of duplication free sequence
#define Num_Ops 1000 // Maximum of the number of operations
#define Max_num 20000 // Maximum of the number of total regions in the dot plot 
#define Max_list 4000 // Maximum of the number of total regions in a list
#define Max_base 500000 // Maximum bases in the dot plot
#define NUM_GAPS 10000 // Maximum of the number of gaps in a dot plot

#define MIN_LEN 500
#define ERR_SM_TH 101
#define ERR_LG_TH 202
#define ERR_TH 305
// #define D_OP_TH 55
#define D_OP_TH 65

#define SMALLEST_TH 3
#define UNIT_THRESHOLD 10
// #define THRESHOLD 50
#define THRESHOLD 5
#define LOOSEN_T 2000
#define E_DIS_TH 2000
#define DIS_THRESHOLD 1500 // Threshold to simplify the dot plot by merging
#define MDIS_THRESHOLD 1000000
#define CLOSE_TH 5000
#define T_DIS_TH 4000 // Threshold for checking the transitivity
#define NEW_DIS_TH 3000 // Threshold to check the alignment which look like almost self alignment
#define NEW_CHECK_TH 50000
#define O_TH 50 // Threshold to be used for checking overlap
#define TD_TH 80 // Threshold to be used for checking overlap
#define TIGHT_TH 1500 // Threshold to use for tight_subset()
#define T_OP_TH 5 // Threshold to check overlapping in deletion
#define RANGE_TH 1000 // Threshold of the range of alignments
// #define RP_BD 220 // Threshold of boundary of a gap and a repeat
#define RP_BD 302 // Threshold of boundary of a gap and a repeat
#define LEN_TH 2000
#define SIZE_OF_ONE_UNIT 5
#define BASIS 3

#define LG_TH 3000 // the standard to discriminate a large alignment and a small alignment 

#define S_RP_BD 100 // Threshold of boundary of a gap and a repeat
#define DEL_TH 0 // Threshold in the decision of the occrrence of deletion
// #define DEL_TH 10 // Threshold in the decision of the occrrence of deletion
// #define DEL_TH 80 // Threshold in the decision of the occrrence of deletion
// #define DEL_TH 20
#define MAX_DEL 5000 // Threshold in the decision of the occrrence of deletion
#define DEL_DIF_TH 100 // Threshold to check whether gaps are from the same deletion or not
#define DEL_ADJ_TH 25 // Threshold to adjust the point of deletion

#define DIF_TH 2 // the threshold to decide the difference between two identities
#define C_TH 10
#define M_THRESHOLD 0
#define M_TH 250 // Threshold in the decision of the occrrence of insertion
#define L_M_TH 500 // Threshold in the decision of the occrrence of insertion
// #define L_M_TH 5000 // Threshold in the decision of the occrrence of insertion
#define T_RATE 0.90
#define RATE_DIF 0.1

#define M_VALUE 800
// #define EFFECTIVE_VALUE 30
// #define EFFECTIVE_VALUE 12
#define EFFECTIVE_VALUE 0
#define NUM_LOOP 10

#define SELF_ID -1

#define G_MODE 0
#define D_MODE 1
#define C_MODE 2 // the combined alignment of two self-alignments and an inter-species alignment 
#define ED_MODE 3 // the coordinates of the combined alignment do not change
#define RD_MODE 4

#define S_FACTOR 1 // the scaling factor
#define NUM_CHARS_LINE 50 // the number of characters per a line

#define LARGE_OVERLAP 100

#define LEFT_SIDE -1
#define TIE 0
#define RIGHT_SIDE 1

#define UNDEF 0
#define ORF 1
#define MATCH 2
#define PARTIAL 3
#define REDUN 4 // appear redundantly

struct I {
	int lower;
	int upper;
};

struct pair {
	int fl; // flag to indicate a pair or a leaf for a left child
	int fr; // flag for a right child
	int lg;
	int rg;
};

struct i_node {
	int gid; // gene id
	int pid; // pair id
};

struct g_list{
	int gid; // a gene id
	int sid; // a species id
	char gname[LEN_LONG_NAME];
	char sname[LEN_NAME];
	char chrname[LEN_LONG_NAME];
	char info[LEN_LONG_NAME];
  char strand; // if the strand is '+', the direction is non-reverse
	          // otherwise, the strand is reverse
	int txStart, txEnd; // the location of a gene
	int cdsStart, cdsEnd; // the location of coding region
	int exonCount;
	int ortho_id; // ortho_id or number of hits
	int type; // 0: undefined, 1: orfs, 2: complete-match, 3: partial-match
//	char *exonReg;
};

struct ops_list{
	char sign; // +/-: dup, d/l: del, c/v: con, s: spe, i: inv
	int id; // id of the starting alignment in the initial dot plots
	int srcStart, srcEnd; // the boundary of the source region
	int dstStart, dstEnd; // the boundary of the target region
	int src_b, src_e; // the ancestral boundary of the source region
	int dst_b, dst_e; // the ancestral boundary of the target region
	float pid; // the percentage identity
	int sp_id; // the species id
};

struct p_tree{
  struct p_tree *left, *right;
	struct p_tree *parent;
	struct I reg;
  char *name; // a gene name or species name
	double b_len; // a branch length
  int d_mode; // if this region is an internal node, d_mode is one of three modes - DP, CV, and SP 
	int od;
  int sp_code; // a code number of self-alignment: species id
  int gid; // a gene identifier 
	int nid; // a node identifier
  int val;
	int *ch_sp; // the list of children species
	int num_csp; // the number of child species
	int depth; // the depth of each node
	bool visited; 
};

struct n_pair{
	char name1[LEN_NAME];	
	char name2[LEN_NAME];	
	char name3[LEN_NAME];	
	int id;
	int len;
};


struct slist{
	int id;
	int val;
	int val_red; // the number of reduced local alignments
	int sp_state;
	int add_sp_state;
	bool is_x;
};

struct ID_List{ // the data structure of the suspend list
	bool is_x; // if is_x is true, x is inserted
	           // otherwise, y is inserted
	int m_id; // id of the inserted region
	int left_id;
	int right_id;
	bool f_is_x; // if f_is_x is true, x of left_id and right_id is involved
	bool is_t_ins; // if is_t_ins is true, the inserted region is tandem
};

struct PR_List{
	int n_match;
	int pid;
	int origin; // if origin is 0, x is duplicated
	            // if origin is 1, y is duplicated
							// if origin is 2, both of x and y are possible as duplicated
							// if origin is 3, it can not be chosen for removal match
};

struct b_list{
  int b1, e1, len1;
  int b2, e2, len2;
  char strand;
  int pid;
};

struct r_list{
	int start;
	int end;
	int len;
	float d_rate;
	char rn[10];
};

struct gap_list{
	int gid;
	int type;
	/*
		0: overlap vs. overlap
		1: overlap vs. no overlap
		2: no overlap vs. overlap
		3: no overlap vs. no overlap
		COPY_INTO_OWN: copy into the own alignment
		11, 21: Tandem overlap vs. no overlap
		12, 22: Tandem no overlap vs. overlap
	*/
	int id1;
	int id2;
	int x1, x2; // matched pair - classifier
	int y1, y2; // gap interval
	int offset;
};

struct DotList{
  struct I x, y;
  int sign; // if the sign is 0, the orientation is the same
	          // else if the sign is 1, the orientation of two regions is reverse
						// else if the sign is 10 or 11, the alignment is temporally hiden
						// else if the sign is 20 or 21, the alignment is involved in merging
						// else if the sign is 2, it would be deleted
	int init_sign; // the initial sign
	struct I m_x, m_y; // the modified region info
	int fid; // ith alignment of the initial dot plot
	int index; // the index of the initial dot plot
	int identity;
	int m_pid; // the original percentage identity
	int l_id; // a linked local alignment split by copy into the own region
            // default is -1
	int l_pid;
	int lock; // a parameter of the locking status, default is -1
	          // if the lock is LEFT_LOCK, a left end can not be merged
						// else if the lock is RIGHT_LOCK, a right end can not be merged
						// else if the lock is BOTH_LOCK, the alignment can not be merged
	int xl_diff, yl_diff, xr_diff, yr_diff;
	int c_id, m_id; // if a local alignment is a chained alignment, this indicates c_id in the chain_list 
						// m_id is another alignment in the chained one 
						// otherwise, -1 by default
	int pair_self; // PAIR or SELF
	int sp_id; // a code of each inter-species and intra-species alignment
	int left_diff, right_diff; // the length of parts cut off while converting the intervals for handling tandem duplication
	int rp1_id, rp2_id; // used for the chaining program and for the tandem dup detection in the duplication inferrence program
	int xl_offset, yl_offset, xr_offset, yr_offset;
};

struct chain_list{
	int c_id; // the identifier in the alignments
	int num_algs; // the number of local alignments in the chain
	int fid[100]; // the list of fid of the initial dot plot in the chain
};

struct BP{
	int pt;
	int reuse_s;
	int reuse_t;
};

struct cv_list{
	int oid; // original id
	int fid; // ith alignment of the initial dot plot
	int s1, s2, t1, t2; // coordinates of the source and target region for conversion
	int a1, a2, b1, b2; // coordinates of the alignment involved in the conversion ; the interval is closed
	int c1, c2, d1, d2;
	int len1, len2;
	char ori, ori1, ori2; // the orientation of the alignment
	float pid;
	char pval[50];
	char name1[50], name2[50], pr_name[50], name3[50];
	char ortho1[50], ortho2[50];
	int sp_id; // 0 for the first species, 1 for the second species
	int dir;
	int bid1, bid2[50]; // branch id
	char bid[500];
	int num_bid2;
	int status;
};

struct exons_list{
	struct I reg;	
	int fid;
	char name[LEN_LONG_NAME];
	char chrname[LEN_LONG_NAME];
};

struct blast_list{
	struct I reg;	
	float pid1, pid2;
	char name1[LEN_LONG_NAME], name2[LEN_LONG_NAME];
	char info1[LEN_LONG_NAME], info2[LEN_LONG_NAME];
};

struct short_alist{ // a reduced structure for a local alignment
	int id; // the id in the original initial dot plot
	struct I x; // a genomic region of species 1
	struct I y;  // a genomic region of species 2 
	int val;
};

struct sp_list{
  char name[LEN_NAME];
  int id;
};

struct scaffold{
	char name[LEN_NAME];
	int b, e;
	int len;
};

void numtostr(int c, char *str);
void output_ops(int num_ops, struct ops_list *ops, int len, float sc);
void assign_dot_list(struct DotList *old_algn, int num_pair, struct DotList *new_algn);
int index_current_sp( int code, int *list1, int *list2);
bool is_current_sp(int code, int *list1, int *list2, int num_sp);

#endif /* MAIN_H */
