// interval: [a,b]
#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"
#include "regions.h"

#define INTERVAL_LEN 100
#define LEN_RATIO 0.2
#define LEN_FRACTION 0.8
#define HOMOLOGOUS_MATCH_PID 75

int debug_mode;
char S[BIG], T[BIG];

bool is_in_orfs_list(char *name, struct sp_list *orfs, int num_orfs);
void mark_in_list(char *name, struct sp_list *orfs, int id, int num_orfs);
int find_id_in_sp_list( char *name, struct sp_list *strains, int num_strains); 

int main(int argc, char **argv)
{
	struct DotList *algns;
	struct DotList *homologs;
	int num_homologs = 0;
	int *num_algns;
	int *size1, *size2;
	FILE *f;
	int i = 0, j = 0;
	int count = 0;
	char species[MAX_NAME], species2[MAX_NAME];
	bool is_end = false;
	int len_cutoff = 0;
	int pid_cutoff = 0;
	int min_len = 0;
	int len1 = 0, len2 = 0;
	struct sp_list *orfs;
	int num_orfs = 0;
	struct sp_list *strains;
	int num_strains = 0;
	int num_clusters = 0;
	int **presence_map;
	int *strain_count;
	int cur_id = -1;
	int temp_count = 0;
	bool is_paralogs = false;

	pid_cutoff = HOMOLOGOUS_MATCH_PID;
	debug_mode = FALSE;
	len_cutoff = INTERVAL_LEN;
	if( argc == 5 ) {
		debug_mode = TRUE;
		len_cutoff = atoi(argv[2]);
		pid_cutoff = atoi(argv[3]);
	}
	else if( argc == 4 ) {
		len_cutoff = atoi(argv[2]);
		pid_cutoff = atoi(argv[3]);
	}
	else if( argc == 3 ) {
		if( strcmp(argv[2], "PARALOGS") == 0 ) {
			is_paralogs = true;
		}
		else {
			pid_cutoff = atoi(argv[2]);
		}
	}
	else if( argc != 2 ) {
		fatal("args: maf\n");
	}

	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	strcpy(S, "");
	strcpy(T, "");

	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[1]);    
	}
	else {
		while((is_end == false) && fgets(S, BIG, f)) {
			if( S[0] == '#' ) {
				while( (is_end == false) && (S[0] == '#') ) {
					if( fgets(S, BIG, f) == NULL ) is_end = true;
					if( strncmp(S, "##maf", 5) == 0 ) {
						count = 0;
					}
				}
				count = 0;	
			}

 		 	if( S[0] == 'a' ) {
				count++;
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", species) != 1) || (sscanf(T, "%*s %s %*s", species2) != 1)) {}
			}
		}
	}
	fclose(f);

	j = 0;
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * (sizeof(struct DotList)) );
		homologs = (struct DotList *) ckalloc(count * (sizeof(struct DotList)) );
		orfs = (struct sp_list *) ckalloc(count * (sizeof(struct sp_list)) );
		strains = (struct sp_list *) ckalloc(count * (sizeof(struct sp_list)) );
 		read_maf(argv[1], G_MODE, algns, num_algns, size1, size2);

		for( i = 0; i < (*num_algns); i++ ) {
			if( strcmp(algns[i].name1, algns[i].name2) == 0 ) {
				if( sscanf(algns[i].name1, "%[^.].%s", species, species2) != 2 ) {
					fatalf("wrong ORF name in %s\n", algns[i].name1);								
				}
				else {
					if( is_in_orfs_list(species, strains, num_strains) == false ) {
						strcpy(strains[num_strains].name, species);
						strains[num_strains].id = -1;
						num_strains++;
					}
				}

				if( sscanf(algns[i].name2, "%[^.].%s", species, species2) != 2 ) {
					fatalf("wrong ORF name in %s\n", algns[i].name2);								
				}
				else {
					if( is_in_orfs_list(species, strains, num_strains) == false ) {
						strcpy(strains[num_strains].name, species);
						strains[num_strains].id = -1;
						num_strains++;
					}
				}

				if( is_in_orfs_list(algns[i].name1, orfs, num_orfs) == false ) {
					strcpy(orfs[num_orfs].name, algns[i].name1);
					orfs[num_orfs].id = -1;
					num_orfs++;
				}

				if( is_in_orfs_list(algns[i].name2, orfs, num_orfs) == false ) {
					strcpy(orfs[num_orfs].name, algns[i].name2);
					orfs[num_orfs].id = -1;
					num_orfs++;
				}
			}
			else if( (is_paralogs == false) && (((float)width(algns[i].x) / (float)algns[i].len1) < LEN_FRACTION)){} 
			else if( (is_paralogs == false) && (((float)width(algns[i].y) / (float)algns[i].len2) < LEN_FRACTION)){} 
			else if( width(algns[i].x) > len_cutoff ) {
				if( algns[i].identity >= pid_cutoff ) {
					len1 = width(algns[i].x);
					len2 = width(algns[i].y);
					if( len1 > len2 ) min_len = len2;
					else min_len = len1;

//					printf("%s %d %d\n", algns[i].name1, algns[i].x.lower, algns[i].x.upper);
					if( (float)((float)abs(len1-len2)/(float)min_len) <= LEN_RATIO ) {
						assign_algn(homologs, j, algns[i]);
						
						if( sscanf(algns[i].name1, "%[^.].%s", species, species2) != 2 ) {
							fatalf("wrong ORF name in %s\n", algns[i].name1);								
						}
						else {
							if( is_in_orfs_list(species, strains, num_strains) == false ) {
								strcpy(strains[num_strains].name, species);
								strains[num_strains].id = -1;
								num_strains++;
							}
						}

						if( sscanf(algns[i].name2, "%[^.].%s", species, species2) != 2 ) {
							fatalf("wrong ORF name in %s\n", algns[i].name2);								
						}
						else {
							if( is_in_orfs_list(species, strains, num_strains) == false ) {
								strcpy(strains[num_strains].name, species);
								strains[num_strains].id = -1;
								num_strains++;
							}
						}

						if( is_in_orfs_list(algns[i].name1, orfs, num_orfs) == false ) {
							strcpy(orfs[num_orfs].name, algns[i].name1);
							orfs[num_orfs].id = -1;
							num_orfs++;
						}
	
						if( is_in_orfs_list(algns[i].name2, orfs, num_orfs) == false ) {
							strcpy(orfs[num_orfs].name, algns[i].name2);
							orfs[num_orfs].id = -1;
							num_orfs++;
						}
	
						j++;
					}
				}
			}
		}
		num_homologs = j;

		num_clusters = 0;
		for( i = 0; i < num_orfs; i++ ) {
			temp_count = 0;
			if( orfs[i].id == -1 ) {
				orfs[i].id = num_clusters;
				for( j = 0; j < num_homologs; j++ ) {
					if( strcmp(orfs[i].name, homologs[j].name1) == 0 ) {
						mark_in_list(homologs[j].name2, orfs, num_clusters, num_orfs);	
						temp_count++;
					}				

					if( strcmp(orfs[i].name, homologs[j].name2) == 0 ) {
						mark_in_list(homologs[j].name1, orfs, num_clusters, num_orfs);	
						temp_count++;
					}
				}
				if( temp_count > 0 ) num_clusters++;
				else {
					orfs[i].id = -1;
				}
			}
			else {
				for( j = 0; j < num_homologs; j++ ) {
					if( strcmp(orfs[i].name, homologs[j].name1) == 0 ) {
						mark_in_list(homologs[j].name2, orfs, orfs[i].id, num_orfs);	
					}				

					if( strcmp(orfs[i].name, homologs[j].name2) == 0 ) {
						mark_in_list(homologs[j].name1, orfs, orfs[i].id, num_orfs);	
					}
				}
			}
		}

		if( ( num_clusters > 0 ) && ( num_strains > 0 ) ) {
			presence_map = (int **) ckalloc(num_clusters * sizeof(int *));
			strain_count = (int *) ckalloc(num_clusters * sizeof(int));
			for( i = 0; i < num_clusters; i++ ) {
				strain_count[i] = 0;
				presence_map[i] = (int *) ckalloc(num_strains * sizeof(int));
				for( j = 0; j < num_strains; j++ ) presence_map[i][j] = 0;
			}
		}

		temp_count = 0;
		for( i = 0; i < num_clusters; i++ ) {
			if( is_paralogs == false ) {
				printf("ORF%d: ", i+1);
			}
			for( j = 0; j < num_orfs; j++ ) {
				if( orfs[j].id == i ) {
					if( is_paralogs == false ) {
						printf("%s ", orfs[j].name);
					}
					temp_count++;
					if( sscanf(orfs[j].name, "%[^.].%s", species, species2) != 2 ) {
					}
					else {
						cur_id = -1;
						cur_id = find_id_in_sp_list( species, strains, num_strains);
						if( cur_id == -1 ) {
							fatalf("%s not in strain list", species);
						}
						else {
							presence_map[i][cur_id] = 1;
						}
					}
				}
			}
			if( is_paralogs == false ) printf("\n");
		}

//		if( ( is_paralogs == false ) || ( (is_paralogs == true) && (num_orfs >= 2) ) ) { 
		j = num_clusters + 1;
		for( i = 0; i < num_orfs; i++ ) {
			if( orfs[i].id == -1 ) {
				if( is_paralogs == false ) { 
						printf("ORF%d: %s\n", j, orfs[i].name);
				}
				temp_count++;
				j++;
			}
		}

		if( is_paralogs == false ) {
			printf("\nTotal %d ORFs\n", temp_count);
			printf("\n");
		}
		for( i = 0; i < num_clusters; i++ ) {
			temp_count = 0;
			for( j = 0; j < num_strains; j++ ) {
				if( presence_map[i][j] == 1 ) temp_count++;
			}
			strain_count[i] = temp_count;	
		}
	
		if( is_paralogs == false ) {
			for( j = 0; j < num_strains; j++ ) printf("%s\t", strains[j].name);
			printf("\n");
			for( i = 0; i < num_clusters; i++ ) {
				if( strain_count[i] > 1 ) {
					for( j = 0; j < num_strains; j++ ) printf("%d\t", presence_map[i][j]);
					printf("\n");
				}
			}
		}

		if( is_paralogs == true ) {
			temp_count = 1;
			for( i = 0; i < num_clusters; i++ ) {
				if( strain_count[i] > 1 ) {
					printf("ORF%d: ", temp_count);
					temp_count++;
					for( j = 0; j < num_orfs; j++ ) {
						if( orfs[j].id == i ) {
							printf("%s ", orfs[j].name);
						}
					}
					printf("\n");
				}
			}
		}
	}	

	free(num_algns);
	if( count > 0 ) {
		if( (num_clusters > 0) && (num_strains > 0) ) {
			for(i = 0; i < num_clusters; i++) {
				free(presence_map[i]);
			}
			free(presence_map);
			free(strain_count);
		}
		free(algns);
		free(homologs);
		free(orfs);
		free(strains);
	}
	free(size1);
	free(size2);
	return EXIT_SUCCESS;
}

int find_id_in_sp_list( char *name, struct sp_list *strains, int num_strains ) 
{
	bool is_in = false;
	int res = -1;
	int i = 0;

	while( (is_in == false) && (i < num_strains) ) {
		if( strcmp(name, strains[i].name) == 0 ) {
			is_in = true;
			res = i;
		}
		i++;
	}

	return(res);
}

bool is_in_orfs_list(char *name, struct sp_list *orfs, int num_orfs)
{
	bool res = false;
	int i = 0;

	while( (res == false) && (i < num_orfs) ) {
		if( strcmp(name, orfs[i].name) == 0 ) res = true;
		i++;
	}

	return(res);
}

void mark_in_list(char *name, struct sp_list *orfs, int id, int num_orfs)
{
	bool res = false;
	int i = 0;

	while( (res == false) && (i < num_orfs) ) {
		if( strcmp(name, orfs[i].name) == 0 ) {
			if( (orfs[i].id == -1) || (id < orfs[i].id) ) {
				orfs[i].id = id;
			}
			res = true;
		}
		i++;
	}
}
