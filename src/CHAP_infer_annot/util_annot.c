#include "util_annot.h"
#include "util.h"

/*
void initialize_exons_list(struct exons_list *a, int from, int to)
{
  int i = 0;
  for( i = from; i < to; i++ ) {
    a[i].reg = assign_I(0,1);
    strcpy(a[i].chr, "");
    strcpy(a[i].name, "");
  }
}
*/

void get_gene_name(char *head, char *gname)
{
	char next[500];
	int len = 0;
	int i = 0, j = 0;

	strcpy(next, "");
	len = strlen(head);

	while((i < len) && (!isspace(head[i])) && ((head[i] != '=')) && (!isspace(head[i]))) i++;

	if( head[i] == '=' ) {
		i++;
		j = i;
		while( (j < len) && (head[j] != '_') && (head[j] != ';') && (!isspace(head[j])) ) 
		{
			next[j-i] = head[j];
			j++;
		} 
		next[j-i] = '\0';
	}
	strcpy(gname, next);
}
