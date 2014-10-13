#ifndef UTIL_I_H
#define UTIL_I_H
#include "sort_genes.h"

int width(struct I temp);
bool overlap(struct I a, struct I b);
bool subset(struct I a, struct I b);
bool proper_subset(struct I a, struct I b);
bool equal(struct I a, struct I b);
struct I assign_I(int lower, int upper);
bool in(int r, struct I reg);
struct I intersect(struct I a, struct I b);

#endif /* UTIL_I_H */
