#ifndef UTIL_GEN_H
#define UTIL_GEN_H
#include "main.h"
#include "kd_tree.h"

void add_symmetric_points(struct DotList *s_dots, struct DotList *dots, int num);
void sort_list(struct slist *st, struct DotList *dots, int num);
void sort_init_algns(struct slist *st, struct DotList *dots, int num, int mode);
void sort_by_yintercept(struct slist *st, struct DotList *dots, int num);
void quick_sort_dec(struct slist *a, int lo, int hi);
void quick_sort_inc(struct slist *a, int lo, int hi);
void quick_sort_inc_int(struct I *a, int lo, int hi);
void quick_sort_inc_alist(struct short_alist *a, int lo, int hi);
void quick_sort_dec_alist(struct short_alist *a, int lo, int hi);
struct slist assign_slist(struct slist a);
struct short_alist assign_alist(struct short_alist a);
void print_sorted(struct slist *a, struct DotList *dots, int num);
void sort_by_pid(struct slist *st, struct DotList *dots, int num);
void sort_by_width(struct slist *st, struct DotList *dots, int num);
void sort_by_width_y(struct slist *st, struct DotList *dots, int num);
void sort_by_loc_short_alist(struct short_alist *list, int num);
void quick_sort_plist_x(struct perm_pt *a, int lo, int hi);
void quick_sort_plist_y(struct perm_pt *a, int lo, int hi);
struct perm_pt assign_pm_val(struct perm_pt a);
void make_into_one(int num, struct DotList *dots, struct DotList *s_dots, struct DotList *t_dots);
int increase_count(struct DotList *org, struct DotList *dots, int id, int c, bool is_x);
void assign_item(struct DotList *dots, int loc, struct I x, struct I y, int id, int s);
void assign_algn(struct DotList *temp, int loc, struct DotList cur);
bool is_match_algn(struct DotList cur, struct DotList cmp);
int find_match_algn(struct DotList cur, struct DotList *algns, int num_algns, struct DotList *temp, int sp_id);
bool merge_into_one_algn(struct DotList *temp, struct DotList cur);
void mark_match_algn(struct DotList cur, struct DotList *algns, int num_algns, int sp_id, int mode);
int assign_sign(int sign, int mode);
bool is_same_sign(int sign1, int sign2);
void sort_rev_init_algns(struct slist *st, struct DotList *dots, int num, int mode);
int search_range_b(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode);
int search_range_e(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode);
int quick_search_range_b(struct slist *sorted, struct DotList *algns, int i, int j, int query, int mode);
int quick_search_range_e(struct slist *sorted, struct DotList *algns, int i, int j, int query, int mode);


#endif /* UTIL_GEN_H */
