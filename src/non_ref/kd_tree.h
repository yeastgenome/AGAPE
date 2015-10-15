#ifndef KD_TREE_H
#define KD_TREE_H

struct kdnode {
	int bucket;
	int cutdim; // if cutdim is 1, the basis is x axis.
							// if it is 2, the basis is y axis
	int cutval;
	struct kdnode *loson, *hison;
	int lopt;
	int hipt;
};

struct perm_pt {
	int id;
	int x_pt;
	int y_pt;
	int sign;
};

struct kdnode *build_kd(struct perm_pt *p_pts, int l, int u);
void assign_perm(struct perm_pt *p_pts, int num, struct DotList *points, int side);
int findmaxspread(struct perm_pt *p_pts, int l, int u);
int adjust_m(struct perm_pt *p_pts, int l, int u, int m, int cutdim);
int px(struct perm_pt *p_pts, int m, int cutdim);
void free_kd(struct kdnode *t);

#endif /* KD_TREE_H */
