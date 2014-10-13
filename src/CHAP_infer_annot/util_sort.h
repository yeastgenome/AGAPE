#ifndef UTIL_SORT_H
#define UTIL_SORT_H

struct slist{
  int id;
  int val;
};

void quick_sort_dec(struct slist *a, int lo, int hi);
void quick_sort_inc(struct slist *a, int lo, int hi);
struct slist assign_slist(struct slist a);

#endif /* UTIL_SORT_H */
