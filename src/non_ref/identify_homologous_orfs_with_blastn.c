// interval: [a,b]
#include "main.h"
#include "util.h"
#include "input_blastn.h"

#define INTERVAL_LEN 100
#define LEN_RATIO 0.2
#define LEN_FRACTION 0.8
#define HOMOLOGOUS_MATCH_PID 75

int debug_mode;

void mark_in_list(char *name, struct sp_list *orfs, int id, int num_orfs);
int find_id_in_sp_list( char *name, struct sp_list *strains, int num_strains); 

int main(int argc, char **argv)
{
	int num_orfs = 0;
	struct sp_list *orfs;
	struct sp_list **homologs;
	struct sp_list *strains;
	int num_strains = 0;
	int num_homologs = 0;
	int num = 0;
	int *max_hits;
	FILE *f;
	int i = 0, j = 0, k = 0;
	int count = 0;
	char species[MAX_NAME], species2[MAX_NAME];
	int len_cutoff = 0;
	int pid_cutoff = 0;
	int num_clusters = 0;
	int **presence_map;
	int *strain_count;
	int *group_count;
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
		fatal("args: blastn_out\n");
	}

	max_hits = (int *) ckalloc(sizeof(int));
	*max_hits = 1;
	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("cannot find alignment in %s", argv[1]);    
	}
	else {
		count = count_genes_blastn(f, max_hits);
	}

	j = 0;
	if( count > 0 ) {
		homologs = (struct sp_list **) ckalloc(count * (sizeof(struct sp_list *)) );
		for( i = 0; i < count; i++ ) {
			homologs[i] = (struct sp_list *) ckalloc((*max_hits) * (sizeof(struct sp_list)));
			for( j = 0; j < (*max_hits); j++ ) {
				homologs[i][j].id = 0;
				strcpy(homologs[i][j].name, "");
			}
		}
		strains = (struct sp_list *) ckalloc(count * (sizeof(struct sp_list)) );
		orfs = (struct sp_list *) ckalloc(count * (sizeof(struct sp_list)));
		for( i = 0; i < count; i++ ) {
			strains[i].id = -1;
			strcpy(strains[i].name, "");
			orfs[i].id = -1;
			strcpy(orfs[i].name, "");
		}
 		num_homologs = input_genes_blastn(f, homologs, HOMOLOGOUS_MATCH_PID, false);

		num_strains = 0;
		num_orfs = 0;
		for( i = 0; i < num_homologs; i++ ) {
			if( sscanf(homologs[i][0].name, "%[^.].%s", species, species2) != 2 ) {
        fatalf("wrong ORF name in %s\n", homologs[i][0].name);
      }
      else {
        if( is_name_already_in(species, strains, num_strains) == false ) {
          strcpy(strains[num_strains].name, species);
          strains[num_strains].id = -1;
          num_strains++;
        }
      }

      if( is_name_already_in(homologs[i][0].name, orfs, num_orfs) == false ) {
        strcpy(orfs[num_orfs].name, homologs[i][0].name);
	      orfs[num_orfs].id = -1;
        num_orfs++;
      }
		}

		num_clusters = 0;
		for( i = 0; i < num_orfs; i++ ) {
			if( orfs[i].id == -1 ) {
				for( j = 0; j < num_homologs; j++ ) {
					if( strcmp(orfs[i].name, homologs[j][0].name) == 0 ) {
						num = homologs[j][0].id;
						if( num > 1 ) {
							orfs[i].id = num_clusters;
							mark_in_list(homologs[j][0].name, orfs, num_clusters, num_orfs);	
							for( k = 1; k < num; k++ ) {
								mark_in_list(homologs[j][k].name, orfs, num_clusters, num_orfs);	
							}
							num_clusters++;
						}
					}				
				}
			}
			else {
				for( j = 0; j < num_homologs; j++ ) {
					if( strcmp(orfs[i].name, homologs[j][0].name) == 0 ) {
						num = homologs[j][0].id;
						for( k = 1; k < num; k++ ) {
							mark_in_list(homologs[j][k].name, orfs, orfs[i].id, num_orfs);	
						}
					}
				}
			}
		}

		if( ( num_clusters > 0 ) && ( num_strains > 0 ) ) {
			presence_map = (int **) ckalloc(num_clusters * sizeof(int *));
			strain_count = (int *) ckalloc(num_clusters * sizeof(int));
			group_count = (int *) ckalloc(num_clusters * sizeof(int));
			for( i = 0; i < num_clusters; i++ ) {
				strain_count[i] = 0;
				group_count[i] = 0;
				presence_map[i] = (int *) ckalloc(num_strains * sizeof(int));
				for( j = 0; j < num_strains; j++ ) presence_map[i][j] = 0;
			}
		}

		for( i = 0; i < num_clusters; i++ ) {
			temp_count = 0;
			for( j = 0; j < num_orfs; j++ ) {
				if( orfs[j].id == i ) {
					temp_count++;
				}	
			}
			group_count[i] = temp_count;
		}

		temp_count = 0;
		k = 0;
		for( i = 0; i < num_clusters; i++ ) {
			if( is_paralogs == false ) {
				if( group_count[i] > 0 ) {
					k++;
					printf("Group%d: ", k);
				}
			}
			for( j = 0; j < num_orfs; j++ ) {
				if( orfs[j].id == i ) {
					if( is_paralogs == false ) {
						printf("%s ", orfs[j].name);
					}
					temp_count++;
					if( sscanf(orfs[j].name, "%[^.].%s", species, species2) != 2 ) {
						fatalf("wrong name in %s\n", orfs[j].name);
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
			if( is_paralogs == false ) {
				if( group_count[i] > 0 ) printf("\n");
			}
		}

//		if( ( is_paralogs == false ) || ( (is_paralogs == true) && (num_orfs >= 2) ) ) { 
//		j = num_clusters + 1;
		j = k + 1;
		for( i = 0; i < num_orfs; i++ ) {
			if( orfs[i].id == -1 ) {
				if( is_paralogs == false ) { 
						printf("Group%d: %s\n", j, orfs[i].name);
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

	free(max_hits);
	if( count > 0 ) {
		if( (num_clusters > 0) && (num_strains > 0) ) {
			for(i = 0; i < num_clusters; i++) {
				free(presence_map[i]);
			}
			free(presence_map);
			free(strain_count);
			free(group_count);
		}
		free(orfs);
		free(strains);
	}
	fclose(f);
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
