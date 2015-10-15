#ifndef UTIL_I_H
#define UTIL_I_H
#include "main.h"

int width(struct I temp);
bool overlap(struct I a, struct I b);
bool subset(struct I a, struct I b);
bool proper_subset(struct I a, struct I b);
bool equal(struct I a, struct I b);
struct I assign_I(int lower, int upper);
bool in(int r, struct I reg);
struct I intersect(struct I a, struct I b);
bool near_equal(struct I a, struct I b, int th);
bool mostly_overlap(struct I a, struct I b, float cut_ratio);
bool proper_overlap(struct I reg1, struct I reg2);

#endif /* UTIL_I_H */
