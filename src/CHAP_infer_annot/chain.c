#include "util.h"
#include "maf.h"
#include "chain.h"
#include "mz_scores.h"
#include "multi_util.h"

extern int MIN_CHAIN;
extern int OVERLAP_THRESHOLD; 

enum {
  INF=1000000000,	// about as big as they get
  HUGECP=100000,      	// limit on number of HSPs
  MAX_GAP=1500,	        // largest gap between HSPs in a chain
  MIN_LINK=1,		// smallest number HSPs in a chain
  SCALE=1500,		// connection penalty is distance/SCALE
};

#define MIN_BASE_PCT 0.3  // smallest pct of bases in a chain

struct hsp_ali *hsp;

#define px(i,j,G) ((j == 0) ? (G->start[G->perm[i]].b1) : (G->start[G->perm[i]].b2))

enum { BUCKET_SIZE=11 };


typedef struct kdnode Kdnode, *Kdptr;

typedef struct best_predecessor bp_t;

static int connect(struct hsp_ali *g, struct hsp_ali *h) {
	int d1 = h->b1 - g->e1, d2 = h->b2 - g->e2, x;
	  if ( h->b1 < g->b1 || h->b2 < g->b2)
	    return INF;

	x = MAX(d1,d2);
	if (x > MAX_GAP)
		return INF;
	else
		return x/SCALE;
}

static int hsp_cmp_beg1(const void *a, const void *b) {
	return ((struct hsp_ali*)a)->b1 - ((struct hsp_ali*)b)->b1;
}

/* --------------------------- build the K-d tree --------------------------- */

#define SWAP(p,q) \
  do {int t; t = G->perm[p]; G->perm[p] = G->perm[q]; G->perm[q] = t;} \
  /*CONSTCOND*/ while(0)

int partition(int L, int U, int cutdim, const gp_t *const G)
{
	int m, a, b, c, v, i, j;

	if (U - L < 2)
		fatal("partition: cannot happen");

	m = (L+U)/2;
	v = a = px(L, cutdim, G);
	b = px(m, cutdim, G);
	c = px(U, cutdim, G);
	// find the median of three entries; move to the front
	if ((a <= b && b <= c) || (c <= b && b <= a)) {
		SWAP(L,m);
		v = b;
	} else if ((a <= c && c <= b) || (b <= c && c <= a)) {
		SWAP(L,U);
		v = c;
	}

	// move smaller entries to front, larger to rear
	i = L;
	j = U+1;
	while (i < j) {
		// search forward for a large entry
		for (++i; i <= U && px(i,cutdim,G) <= v; ++i)
			;
		// search backward for a small entry
		for (--j; j >= L && px(j,cutdim,G) > v; --j)
			;
		SWAP(i,j);
	}
	SWAP(i,j);	// reverse the last swap
	SWAP(L,j);	// move the "pivot" value to the proper location

	// warning: if we return j = U, build() will recurse forever
	if (j < U)
		return j;
	if (U-L == 2) return U-1;
	return partition(L, U-1, cutdim, G);
}

Kdptr build(int L, int U, int toggle, const gp_t *const G)
{
	Kdptr p;
	int m;

	p = ckalloc(sizeof(Kdnode));
	p->max_score = 0;
	if (U-L+1 <= BUCKET_SIZE) {
		p->bucket = 1;
		p->lopt = L;
		p->hipt = U;
		p->loson = p->hison = NULL;
	} else {
		p->bucket = 0;
		m = partition(L, U, toggle, G);
		p->cutval = px(m, toggle, G);
		p->hipt = m;
		p->loson = build(L, m, 1-toggle, G);
		p->hison = build(m+1, U, 1-toggle, G);
	}
	return p;
}

int overlap_limit(struct hsp_ali* s, struct hsp_ali* t) {
  int over_beg, over_end, size1, size2;
  
  size1 = s->e1 - s->b1 + 1;
  size2 = t->e1 - t->b1 + 1;
  over_beg = s->b1 > t->b1 ? s->b1 : t->b1;
  over_end = s->e1 < t->e1 ? s->e1 : t->e1;
  if ( (double)(over_end - over_beg+1)/size1 > (double)OVERLAP_THRESHOLD/100)
    return 1;
  if ( (double)(over_end - over_beg+1)/size2 > (double)OVERLAP_THRESHOLD/100)
    return 1;
  size1 = s->e2 - s->b2 + 1;
  size2 = t->e2 - t->b2 + 1;
  over_beg = s->b2 > t->b2 ? s->b2 : t->b2;
  over_end = s->e2 < t->e2 ? s->e2 : t->e2;
  if ( (double)(over_end - over_beg+1)/size1 > (double)OVERLAP_THRESHOLD/100)
    return 1;
  if ( (double)(over_end - over_beg+1)/size2 > (double)OVERLAP_THRESHOLD/100)
    return 1;
  return 0;
}

/* -------------------- find a best predecessor of an HSP ------------------- */
bp_t best_pred(Kdptr Q, int lowbd_dist, int toggle,
			bp_t bp, const gp_t *const G)
{
	int i, d;
	big_score_t contrib;
	
	if (Q == NULL)
		fatal("best_pred: Q is NULL");

	if (Q->max_score - lowbd_dist/SCALE <= bp.contrib)
		return bp;

	if (Q->bucket) {
		for (i = Q->lopt; i <= Q->hipt; ++i) {
			int j = G->perm[i];
			struct hsp_ali *s = G->start+j;
			if (s->b1 >= G->X || s->b2 >= G->Y)
			  continue;
			if ( overlap_limit(s, G->query) != 0)
			  continue;
			contrib =
			   s->total_score - (big_score_t)connect(s, G->query);
			if (contrib >= bp.contrib) {
				bp.contrib = contrib;
				bp.num = j;
			}
		}
	} else {
		d = (toggle ? G->Y : G->X) - Q->cutval;
		if (d >= 0)
	      		bp = best_pred(Q->hison, lowbd_dist, 1-toggle, bp, G);
		if (d <= MAX_GAP)
			bp = best_pred(Q->loson, MAX(lowbd_dist, d),
					1-toggle, bp, G);
	}
	return bp;
}

void update_max_score(Kdptr Q, big_score_t score, int pos)
{
	if (Q != NULL) {
		Q->max_score = MAX(Q->max_score, score);
		if (Q->hipt >= pos)
			update_max_score(Q->loson, score, pos);
		else
			update_max_score(Q->hison, score, pos);
	}
}
bp_t null_bp(bp_t bp) {
	bp.num = -1;
	bp.contrib = 0;
	return bp;
}

int make_chain(Kdptr root, gp_t G, int n)
{
	big_score_t best, ts;
	int best_end, i;
	struct hsp_ali *h;
	bp_t bp;

//fprintf(stderr, "make_chain with nhsp = %d\n", n);

	best_end = -1;
	best = 0;
	for (i = 0; i < n; i++) {
		h = G.start + i;
		if (h->used) {
			h->total_score = 0;
			h->chain = -1;
			continue;
		}
		G.query = G.start+i;
		G.X = G.query->b1;
		G.Y = G.query->b2;
		bp = best_pred(root, 0, 1, null_bp(bp), &G);
		ts = h->total_score = (big_score_t)(h->score) + bp.contrib;
		if (ts > best) {
			best = ts;
			best_end = i;
		}
		h->chain = bp.num;
		if ( bp.num != -1) 
		  (G.start + bp.num)->succ = i;
		if (h->chain >= i) 
		  fatal("chain out of order");
		update_max_score(root, ts, G.rev_perm[i]);
	}
	return 0;
}

struct chain* extract_all_chain(gp_t G, int n) {
  struct chain *chn, *head=NULL;
  struct hsp_ali *h1, *h2, *h;
  int i, j, best;
  
  for (i=0; i<n; i++) {
    h1 = G.start + i;
    if ( h1->used == 1)
      continue;
    if ( h1->chain != -1 && (G.start+h1->chain)->used != 1 )
      continue;
    for (j=i; j<n; ) {
      h2 = G.start + j;
      if ( h2->succ == -1 || (G.start+h2->succ)->used == 1 )
	break;
      j = h2->succ;
    }
    chn = (struct chain*)malloc(sizeof(struct chain));
    chn->score = 0;
    chn->b1 = hsp[i].b1;  chn->b2 = hsp[i].b2;
    chn->e1 = hsp[j].e1;  chn->e2 = hsp[j].e2;
    chn->len1 = chn->len2 = chn->keep_len = 0;
    chn->used_array_1 = chn->used_array_2 = NULL;
    chn->segments = NULL;
    for (best = j; best>=0 && hsp[best].used==0; best=hsp[best].chain ) {
      hsp[best].next = chn->segments;
      chn->segments = &(hsp[best]);
      hsp[best].used = 1;
      chn->len1 += hsp[best].ali->components->size;
      chn->len2 += hsp[best].ali->components->next->size;
      chn->score += hsp[best].score;
    }
    for (h=chn->segments; h!=NULL; h=h->next)
      h->ali->chain_len = chn->len1;
    chn->next = head;
    chn->contigA = copy_string(hsp[i].ali->components->src);
    chn->contigB = copy_string(hsp[i].ali->components->next->src);
    head = chn;
  }
  return head;
}

/*
struct chain* extract_chain(Kdptr root, gp_t G, int n) {
  big_score_t best;
  int best_end, i,j, nlink, chain_len, bases, last;
  struct hsp_ali *h;
  struct chain* ptrChain;

        for (i=0, best=0, best_end=-1; i<n;i++) {
	  h = G.start + i;
	  //if ( h->used )
	  //continue;
	  if ( !h->used && h->total_score > best) {
	    best = h->total_score;
	    best_end = i;
	  }
	}
	if ( best_end < 0)
	  return NULL;
	
	for (j = best_end, nlink=0, bases=0; j != -1; j = hsp[j].chain) {
	  //	  hsp[j].used = 1;
	  ++nlink;
	  bases += (hsp[j].e1-hsp[j].b1+1);
	  if (hsp[j].chain == -1 || hsp[hsp[j].chain].used == 1) {
	    chain_len = hsp[best_end].e1 - hsp[j].b1 + 1;
	    last = j;
	    break;
	  }
	}
	ptrChain = (struct chain*)malloc(sizeof(struct chain));
	ptrChain->score = 0;
	ptrChain->b1 = hsp[last].b1;
	//if ( ptrChain->b1 == 9555)
	//printf("hello");
	ptrChain->b2 = hsp[last].b2;
	ptrChain->e1 = hsp[best_end].e1;
	ptrChain->e2 = hsp[best_end].e2;
	//ptrChain->keep_len = ptrChain->e1-ptrChain->b1+1 + ptrChain->e2-ptrChain->b2+1;
	ptrChain->len1 = ptrChain->len2 = ptrChain->keep_len = 0;
	ptrChain->used_array_1 = ptrChain->used_array_2 = NULL;
	ptrChain->segments = NULL;
	ptrChain->next = NULL;
	while ( best_end >= 0 && hsp[best_end].used==0) {
	  hsp[best_end].next = ptrChain->segments;
	  ptrChain->segments = &(hsp[best_end]);
	  hsp[best_end].used = 1;
	  ptrChain->len1 += hsp[best_end].ali->components->size;
	  ptrChain->len2 += hsp[best_end].ali->components->next->size;
	  ptrChain->score += hsp[best_end].score;
	  best_end = hsp[best_end].chain;
	}
	for (h=ptrChain->segments; h!=NULL; h=h->next)
	  h->ali->chain_len = ptrChain->len1;

	return ptrChain;
}
*/

void reset(Kdptr Q) {
	if (Q == NULL)
		return;
	reset(Q->loson);
	reset(Q->hison);
	Q->max_score = 0;
}

void freeKd(Kdptr Q) {
	if (Q == NULL)
		return;
	freeKd(Q->loson);
	freeKd(Q->hison);
	free(Q);
}

struct chain* all_chains(int n) {
	int i;
	Kdptr root;
	gp_t G;
	struct chain* chain_root=NULL;

	G.start = hsp;
	G.perm = ckalloc(sizeof(int)*n);
	G.rev_perm = ckalloc(sizeof(int)*n);

	for (i = 0; i < n; ++i)
		G.perm[i] = i;
//fprintf(stderr, "start buiding the tree ... ");
	root = build(0, n-1, 1, &G);
//fprintf(stderr, "done\n");
	for (i = 0; i < n; ++i) 
	     G.rev_perm[G.perm[i]] = i;
	make_chain(root, G, n);
	
	/*
	while ( (chain_ptr=extract_chain(root, G, n))!=NULL) {
	  chain_ptr->next = chain_root;
	  chain_root = chain_ptr;
	}
	*/
	chain_root = extract_all_chain(G, n);

	freeKd(root);
	root = NULL;

	free(G.perm);
	free(G.rev_perm);

	return chain_root;
}

void print_chain(struct chain* ptr) {
  struct hsp_ali* phsp;
  printf("chain: %d-%d\n", ptr->b1, ptr->e1);
  for (phsp=ptr->segments; phsp!=NULL; phsp=phsp->next) 
    mafWrite(stdout, phsp->ali);
}
struct chain* process_chain_strand(struct mafAli** pphead, char strand) {
  char contigA[2000], contigB[2000];
  struct mafAli *ali, *wk_list, *cp_list, *head=*pphead, *backup=NULL;
  struct hsp_ali *h;
  int i, nhsp;
  struct chain* chain_root, *ret_root=NULL, *chain_ptr;

  while (head != NULL) {
    cp_list = NULL;
    strcpy(contigA, head->components->src);
    strcpy(contigB, head->components->next->src);
    wk_list = head;
    head = head->next;
    wk_list->next = NULL;
    while (head != NULL) {
      ali = head;
      head = head->next;
      if ( strcmp(ali->components->src, contigA) == 0 && strcmp(ali->components->next->src, contigB)==0 ) {
	ali->next = wk_list;
	wk_list = ali;
      }
      else {
	ali->next = cp_list;
	cp_list = ali;
      }
    }
    head = cp_list;

    for (nhsp=0, ali=wk_list; ali!=NULL; ali=ali->next)
      if ( ali->components->next->strand == strand)
	nhsp++;
    hsp = (struct hsp_ali*)malloc(nhsp*sizeof(struct hsp_ali));
    for (i=0; i<nhsp; i++)
      hsp[i].next = NULL;
  
    i = 0;
    for (ali=wk_list; ali!=NULL;ali=ali->next )
      if ( ali->components->next->strand == strand) {
	h=&(hsp[i]);
	h->b1 = ali->components->start;
	h->b2 = ali->components->next->start;
	h->e1 = ali->components->start + ali->components->size - 1;
	h->e2 = h->b2 + ali->components->next->size - 1;
	h->score = ali->score;
	h->used = 0;
	h->chain = h->succ = -1;
	h->ali = ali;
	++i;
      }
    
    qsort((char*)hsp, nhsp, sizeof(struct hsp_ali), hsp_cmp_beg1);    
    chain_root = all_chains(nhsp);

    for (ali=wk_list; ali->next != NULL; ali=ali->next) 
      ; 
    ali->next = backup; 
    backup = wk_list; 

    if ( chain_root == NULL)
      continue;
    for (chain_ptr=chain_root; chain_ptr->next != NULL; chain_ptr=chain_ptr->next)
      ;
    chain_ptr->next = ret_root;
    ret_root = chain_root;
    
  }

  *pphead = backup;
  return ret_root;
}

struct chain* process_chain(struct mafAli* head) {
  struct chain* pos, *neg, *chain_root, *ptr;

  pos = process_chain_strand(&head, '+');
  neg = process_chain_strand(&head, '-');
  if ( neg == NULL) {
    chain_root = pos;
    return chain_root;
  }
  for (ptr = neg; ptr->next!=NULL; ptr=ptr->next)
    ;
  ptr->next = pos;
  chain_root = neg;
  
  return chain_root;
}

int compar_chain_score(const void* a, const void* b) {
  return (*((struct chain **)b))->score - (*((struct chain **)a))->score;
}

void insert_hali_head_chain(struct chain* chn, struct hsp_ali* hali) {
  hali->next = chn->segments;
  chn->segments = hali;
}

void insert_head_chain(struct chain** head, struct chain* target) {
  target->next = *head;
  *head = target;
}

// get score and len between beg and end
int part_chain_score(struct chain* chn, int beg, int end) {
  struct hsp_ali* hali;
  int cbeg, cend, over_beg, over_end, small, large, score=0;

  small = MAX_INT;
  large = MIN_INT;

  for (hali = chn->segments; hali != NULL; hali=hali->next) {
    if ( hali->b1 > end || hali->e1 < beg )
      continue;
    over_beg = hali->b1 > beg ? hali->b1 : beg;
    over_end = hali->e1 < end ? hali->e1 : end;
    if ( over_beg < small)
      small = over_beg;
    if ( over_end > large)
      large = over_end;

    cbeg = mafPos2Col(hali->ali->components, over_beg, hali->ali->textSize);
    cend = mafPos2Col(hali->ali->components, over_end, hali->ali->textSize);
    
    score += (int)mafScoreRange(hali->ali, cbeg, cend-cbeg+1);
  }
  return score;
}

int top2botChainPos(struct chain* chn, int pos, int after) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret, col;

  for (hali=chn->segments; hali!=NULL; hali=hali->next) {
    comp = hali->ali->components;
    if ( pos >= comp->start && pos <= comp->start + comp->size - 1 )
      break;
  }

  if ( hali != NULL) {
    col = mafPos2Col(comp, pos, hali->ali->textSize); 
    if ( after == 1)
      ret = colPos2Maf_after(comp->next, col);
    else
      ret = colPos2Maf_before(comp->next, col);

    if (ret >= 0)
      return ret;
  }
  return -1;
}

int bot2topChainPos(struct chain* chn, int pos, int after) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret, col;

  for (hali=chn->segments; hali!=NULL; hali=hali->next) {
    comp = hali->ali->components->next;
    if ( pos >= comp->start && pos <= comp->start + comp->size - 1)
      break;
  }
  if ( hali != NULL) {
    col = mafPos2Col(comp, pos, hali->ali->textSize);
    if ( after == 1)
      ret = colPos2Maf_after(hali->ali->components, col);
    else
      ret = colPos2Maf_before(hali->ali->components, col);
    if ( ret >= 0)
      return ret;
  }
  return -1;
}

int top2botChainPos_after(struct chain* chn, int pos) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret;

  ret = top2botChainPos(chn, pos, 1);
  if ( ret >= 0)
    return ret;

  ret = MAX_INT;
  for (hali = chn->segments; hali!=NULL; hali=hali->next) {
    comp = hali->ali->components;
    if ( comp->start < pos )
      continue;
    if ( comp->next->start < ret )
      ret = comp->next->start;
  }
  if ( ret == MAX_INT)
    return -1;
  return ret;
}

int bot2topChainPos_after(struct chain* chn, int pos) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret;

  ret = bot2topChainPos(chn, pos, 1);
  if ( ret >= 0)
    return ret;

  ret = MAX_INT;
  for (hali = chn->segments; hali != NULL; hali = hali->next) {
    comp = hali->ali->components;
    if ( comp->next->start < pos)
      continue;
    if ( comp->start < ret )
      ret = comp->start;
  }
  if ( ret == MAX_INT )
    return -1;
  return ret;
}

int top2botChainPos_before(struct chain* chn, int pos) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret;

  ret = top2botChainPos(chn, pos, 0);
  if ( ret >= 0)
    return ret;

  for ( hali = chn->segments; hali!=NULL; hali=hali->next) {
    comp = hali->ali->components;
    if ( comp->start + comp->size - 1 >= pos )
      continue;
    if ( comp->next->start + comp->next->size - 1 > ret)
      ret = comp->next->start + comp->next->size - 1;
  }
  return ret;
}

int bot2topChainPos_before(struct chain* chn, int pos) {
  struct hsp_ali* hali;
  struct mafComp* comp;
  int ret;

  ret = bot2topChainPos(chn, pos, 0);
  if ( ret >= 0)
    return ret;

  for ( hali = chn->segments; hali!=NULL; hali=hali->next) {
    comp = hali->ali->components;
    if ( comp->next->start + comp->next->size - 1 >= pos )
      continue;
    if ( comp->start + comp->size - 1 > ret)
      ret = comp->start + comp->size - 1;
  }
  return ret;
}

struct chain* mafs2chain(int num, char** strvec) {
  struct chain* chain_root=NULL, *chain_ptr, *chn;
  struct mafFile* maf;
  int i;
  for (i=0; i<num; i++) {
    maf = mafReadAll(strvec[i], 0);
    chn = process_chain(maf->alignments);
    if ( chn==NULL)
      continue;
    for (chain_ptr = chn; chain_ptr->next != NULL; chain_ptr = chain_ptr->next)
      ;
    chain_ptr->next = chain_root;
    chain_root = chn;
  }
  return chain_root;
}

/*
//--------< test entry point >-------------------------
int main(int argc, char **argv) {
  FILE* fp;
  struct chain* chain_root, *chain_ptr;
  struct mafFile* maf;
  struct mafAli* head;

	if (argc > 2)
		fatal("arg = nblastz.out (else stdin)");
	maf = mafReadAll(argv[1], 1);
	head = maf->alignments;
	maf->alignments = NULL;
	mafFileFree(&maf);


	chain_root = process_chain(head);
	
	for (chain_ptr=chain_root; chain_ptr != NULL; chain_ptr=chain_ptr->next)
	  print_chain(chain_ptr);
	
	return 0;
}
*/

