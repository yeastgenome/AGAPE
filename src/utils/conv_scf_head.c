#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"

#define BIG 10000
#define TRUE 1
#define FALSE 2

int is_number(char *str);
int main(int argc, char *argv[])
{
	FILE *f;
	char buf[BIG];
	char next[100];
	char scf_prefix[100];
	char scf_name[100];
	int len = 0;
	int max_scf_num = 0;
	int len_prefix = 0;
	int cur_len = 0, cur_num = 0;
	char scf_len_char[100];
	int count = 0;
	int *scf_len;
	int i = 0;

	strcpy(buf, "");
	strcpy(next, "");
	strcpy(scf_prefix, "");
	strcpy(scf_name, "");
	if( argc != 3 ) {
		printf("conv_scf_head seq_file primary_prefix\n");
		return 1;
	}

	strcpy(scf_prefix, argv[2]);
	len_prefix = strlen(scf_prefix);
	strcpy(scf_len_char, "UNDEF");

	f = fopen(argv[1], "r");

	while(fgets(buf, BIG, f)) {
		if(buf[0] == '>') {
			count++;
		}
	}
	
	fseek(f, 0, SEEK_SET);

	if( count > 0 ) {
		scf_len = (int *) ckalloc(count * sizeof(int));
		for( i = 0; i < count; i++ ) scf_len[i] = 0;
	}

	i = -1;
	while(fgets(buf, BIG, f))
	{
		if(buf[0] == '>') 
		{
			if( i >= 0 ) scf_len[i] = count;	
			count = 0;
			i++;

			if( (sscanf(buf+1, "%s %s %*s", next, scf_len_char) == 2) || (sscanf(buf+1, "%s %*s", next) == 1) ) {
				if( strstr(next, scf_prefix) ) {
					strcpy(scf_name, next+len_prefix);
					cur_len = strlen(next)-len_prefix;
					scf_name[cur_len] = '\0';
				
					cur_num = atoi(scf_name);
					if( cur_num > max_scf_num ) {
						max_scf_num = cur_num;
					}
				}	
			}
			else {
				printf("bad format in %s\n", buf);
				return EXIT_FAILURE;
			}
		}
		else {
			len = strlen(buf);
			if( buf[len-1] == '\n' ) count = count + len-1;
			else count = count + len;
		}
	}
	if( i >= 0 ) scf_len[i] = count;	

	len = 0;
	cur_num = max_scf_num + 1; 
	fseek(f, 0, SEEK_SET);
	i = 0;
	while(fgets(buf, BIG, f))
	{
		if( buf[0] == '>' ) {
			strcpy(scf_len_char, "UNDEF");
			if( (sscanf(buf+1, "%s %s %*s", next, scf_len_char) == 2) || (sscanf(buf+1, "%s %*s", next) == 1) ) {
				if( strstr(next, scf_prefix) ) {
					len = strlen(buf);
					if( buf[len-1] == '\n' ) buf[len-1] = '\0';
					if( scf_len[i] != 0 ) {
						sprintf(scf_name, "%s %d", buf, scf_len[i]);
					}
					if( scf_len[i] != 0 ) {
						printf("%s\n", scf_name);
					}
				}
				else {
					if(is_number(scf_len_char) == TRUE) {
						len = atoi(scf_len_char);
						if( len != scf_len[i] ) {
							fatalf("scaffold length counting error %d-%d\n", len, scf_len[i]);
						}
						sprintf(scf_name, "%s%d %d", scf_prefix, cur_num, len);	
					}
					else sprintf(scf_name, "%s%d %d", scf_prefix, cur_num, scf_len[i]);	

					cur_num++;
					if( scf_len[i] != 0 ) {
						printf(">%s\n", scf_name);
					}
				}	
				len = 0;
			}
			else 
			{
				printf("bad format in %s\n", buf);
				return 1;
			}
			i++;
		}
		else {
			len = strlen(buf);
			if( scf_len[i-1] != 0 ) {
				if( buf[len-1] == '\n' ) printf("%s", buf);
				else printf("%s\n", buf);
			}
		}
	}

	if( count > 0) free(scf_len);
	fclose(f);

	return EXIT_SUCCESS;
}

int is_number(char *str)
{
	int i = 0;
	int res = TRUE;

	while( (str[i] != '\0') && (str[i] != '\n') ) {
		if( isdigit(str[i]) == 0 ) res = FALSE;
		i++;
	}

	return(res);
}
