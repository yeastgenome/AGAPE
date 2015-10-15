#define EXTRACT_GENE_CLUSTER
#ifdef EXTRACT_GENE_CLUSTER

#include "chain.h"

struct pos_pair {
  char* contigA, *contigB;
  int start, end, count, ystart, yend;
  struct pos_pair* next;
};

//struct pos_pair* extract_clusters(struct pos_pair* pp_head, int count);
struct pos_pair* extract_clusters(struct pos_pair* pp_head);
struct pos_pair* get_clusters(struct chain*);

#endif
