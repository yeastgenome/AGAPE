#ifndef REGIONS_H
#define REGIONS_H

#include "main.h"

bool proper_in(int r, struct I reg);
bool proper_overlap(struct I reg1, struct I reg2);
bool f_loose_overlap(struct I reg1, struct I reg2, int th);
bool left_loose_overlap(struct I reg1, struct I reg2, int th);
bool right_loose_overlap(struct I reg1, struct I reg2, int th);
bool loose_overlap(struct I reg1, struct I reg2);
bool fully_subset(struct I reg1, struct I reg2);
bool almost_subset(struct I reg1, struct I reg2);
bool too_loosen_subset(struct I reg1, struct I reg2);
bool f_loosen_subset(struct I reg1, struct I reg2, int th);
bool loosen_subset(struct I reg1, struct I reg2);
bool tight_overlap(struct I reg1, struct I reg2);
bool strict_overlap(struct I reg1, struct I reg2, int threshold);
void init_array(int *array, int num);
void overwrite_dots(int *num, struct DotList *dots);
bool almost_equal(struct I reg1, struct I reg2);
bool strict_almost_equal(struct I reg1, struct I reg2);
bool strict_subset(struct I reg1, struct I reg2);
int compute_distance(struct I x1, struct I y1, struct I x2, struct I y2, int sign);

#endif /* REGIONS_H */
