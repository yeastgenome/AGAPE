/*************************************************
 * Extract gene clusters from pairwise alignment *
 * files.                                        *
 * args: reference's self-alignment and alignment*
 * vs other spieces, suppose all topped by ref.  *
 *************************************************/
#include <stdlib.h>
#include "util.h"
#include "maf.h"
#include "chain.h"
#include "multi_util.h"
#include "extract_gene_cluster.h"

extern int MIN_CLUSTER_CHAIN;
extern int OVERLAP_THRESHOLD;

static int compar_start(const void *a, const void *b) {
  return (*((struct pos_pair **)a))->start -
         (*((struct pos_pair **)b))->start;
}

static int compar_end(const void *a, const void *b) {
  return (*((struct pos_pair **)a))->end -
         (*((struct pos_pair **)b))->end;
}

static struct pos_pair* new_pos_pair(char* contigA,char* contigB, int start, int end, int count, int ystart, int yend)  {
  struct pos_pair* pp;
  pp = (struct pos_pair*)malloc(sizeof(struct pos_pair));
  pp->contigA = copy_string(contigA);
  pp->contigB = copy_string(contigB);
  pp->start = start;   pp->end = end;   pp->count = count;
  pp->ystart = ystart; pp->yend = yend; pp->next = NULL;
  return pp;
}

//struct pos_pair* extract_clusters(struct pos_pair* pp_head, int count) {
struct pos_pair* extract_clusters(struct pos_pair* pp_head) {
  char contig[2000];
  int i = 0, j = 0, k = 0, beg = 0, end = 1, flag_cluster = 0, cluster_beg = 0, cluster_end = 1, num = 0;
  struct pos_pair *pp, **pp_start_array, **pp_end_array, *clusters, *ret_cluster=NULL, *cp_list, *wk_list;
  int overlap_count = 0, nonoverlap_count = 0;

  while ( pp_head != NULL) {
    wk_list = pp_head;
    cp_list = NULL;
    clusters = NULL;
    pp_head = pp_head->next;
    wk_list->next = NULL;
    strcpy(contig, wk_list->contigA);
    num=1;
    while (pp_head != NULL) {
      pp = pp_head;
      pp_head = pp_head->next;
      if ( strcmp(contig, pp->contigA)==0) {
	pp->next = wk_list;
	wk_list = pp;
	num++;
      }
      else {
	pp->next = cp_list;
	cp_list = pp;
      }
    }
    pp_start_array = (struct pos_pair**)malloc((num+1)*sizeof(struct pos_pair*));
    pp_end_array = (struct pos_pair**)malloc((num+1)*sizeof(struct pos_pair*));
    for (i=0, pp=wk_list; i<num; i++, pp=pp->next) 
      pp_start_array[i] = pp_end_array[i] = pp;
      
    for (i=0; i<num; i++)
      pp_start_array[i]->next = NULL;
    pp = new_pos_pair(contig, contig, MAX_INT - 3000, MAX_INT, 0, 0, 0);
    pp_start_array[num] = pp_end_array[num] = pp;
    
    qsort((void*)pp_start_array, num+1, sizeof(struct pos_pair*), compar_start);
    qsort((void*)pp_end_array, num+1, sizeof(struct pos_pair*), compar_end);

    beg = MIN_INT;  end = MAX_INT;
    flag_cluster = 0;  cluster_beg = pp_start_array[0]->start;
    for (i=0; i<=num; i++) {
      if ( pp_start_array[i]->start > end ) {
	nonoverlap_count = overlap_count = 0;
	for (j=i-1; j>=flag_cluster; j--) {
	  for (k=j+1; k<i; k++)
	    if ( overlap(pp_start_array[j]->ystart, pp_start_array[j]->yend, pp_start_array[k]->ystart, pp_start_array[k]->yend) ) {
	      overlap_count++;
	      break;
	    }
	  if ( k==i)
	    nonoverlap_count++;
	}
	for (j=i-1, cluster_end = end; j>= 0; j-- ) {
	  if ( pp_start_array[j]->end < pp_start_array[i]->start ) {
	    if ( pp_start_array[j]->end > cluster_end )
	      cluster_end = pp_start_array[j]->end;
	  }
	  else
	    break;
	}
	if ( j>=0 && i!= num)
	  cluster_end = (cluster_end + pp_start_array[i]->start)/2;
	pp = new_pos_pair(contig, contig, cluster_beg, cluster_end, nonoverlap_count, 0, 0);
	pp->next = clusters;   clusters = pp;

        if ( j>=0 && i != num) 
          cluster_beg = cluster_end + 1; 
        else 
          cluster_beg = pp_start_array[i]->start; 
       
	end = pp_start_array[i]->end; // ???

	for (j=i-1; j>=0; j--) 
	  if (pp_start_array[j]->end > cluster_beg && pp_start_array[j]->end < end )
	    end = pp_start_array[j]->end;

	flag_cluster = i;
      } // if
      beg = pp_start_array[i]->start;
      end = pp_start_array[i]->end < end ? pp_start_array[i]->end : end;
    } // for
    
    for (pp=clusters; pp->next != NULL; pp=pp->next)
      ;
    pp->next = ret_cluster;
    ret_cluster = clusters;
    pp_head = cp_list;
  }

  return ret_cluster;
}

struct pos_pair* get_clusters(struct chain* chain_root) {
  int count=0;
  struct chain *chain_ptr;
  struct pos_pair *pp, *pp_head=NULL, *clusters;

  for (chain_ptr=chain_root; chain_ptr != NULL; chain_ptr = chain_ptr->next) {
    if ( chain_ptr->e1 - chain_ptr->b1 + 1 < MIN_CLUSTER_CHAIN || chain_ptr->e2 - chain_ptr->b2 + 1 < MIN_CLUSTER_CHAIN)
      continue;
    pp = new_pos_pair(chain_ptr->contigA, chain_ptr->contigB, chain_ptr->b1, chain_ptr->e1, 0, chain_ptr->b2, chain_ptr->e2);
    pp->next = pp_head; pp_head = pp; count++;
  }
  
//  clusters = extract_clusters(pp_head, count);
	clusters = extract_clusters(pp_head);
  return clusters;
}

#define EXTRACT_TEST
#ifdef EXTRACT_TEST
//----------------< test entry >----------------------
int main(int argc, char** argv) {
  struct pos_pair *clusters, *pp;
  struct chain* chain_root;

  chain_root = mafs2chain(argc-1, argv+1);
  clusters = get_clusters(chain_root);
  for ( pp = clusters; pp!=NULL; pp=pp->next)
    //    if ( pp->count >= 3 && pp->count <= 20)
			printf("%s %d - %d (%d)\n", pp->contigA, pp->start, pp->end, pp->count);
  return EXIT_SUCCESS;
}
//------------------< end of test >-------------------
#endif
