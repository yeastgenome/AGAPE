#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void numtoroman(char *from, char *to)
{
	int num = 0;
	char symbols[30][50] = { "0", "I", "II", "III", "IV", "V", "VI", "VII", 
		"VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII",
		"XIX", "XX", "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII",
		"XXVIII", "XXIX" };

	if( isdigit(from[0]) ) {
		if( from[0] == '0' ) {
			sscanf(from+1, "%d", &num);
		}
		else {
			sscanf(from, "%d", &num);
		}	
		strcpy(to, "chr");
		strcat(to, symbols[num]);
	}
	else strcpy(to, from); 
}

void get_chr_name(char *head, char *chr)
{
	char res[500], roman[500];
	char next[500];
	int len = 0;
	int i = 0, j = 0;

	strcpy(next, "");
	strcpy(res, "");
	strcpy(roman, "");

	while((head[i] != '_') && (!isspace(head[i]))) i++;

	if( head[i] == '_' ) {
		i++;
		for( j = i; j < (i+3); j++ ) next[j-i] = (char) toupper(head[j]);
//strncpy(head+i, next, 3*sizeof(char));
		next[3] = '\0';

		if( strcmp(next, "CHR") == 0 ) {
			i = i+3;
			strcpy(res, head+i);
			numtoroman(res, roman);
		}
		else {
			strcpy(res, head+i);
			if( strcmp(res, "mito") == 0 ) {
				strcpy(roman, "chrM");
			}
			else strcpy(roman, res);
		} 
	}
	else {
		strcpy(roman, head);	
	}

	len = strlen(roman);
	if( roman[len-1] == '\n' ) {
		roman[len-1] = '\0';
	}	

	strcpy(chr, roman);
}

int main(int argc, char *argv[])
{
	FILE *f;
	char buf[500];
	char next[500];
	char chr[500];
	int len = 0;

	if( argc != 2 ) {
		printf("conv_head seq_file\n");
		return 1;
	}

	f = fopen(argv[1], "r");

	while(fgets(buf, 500, f))
	{
		if(buf[0] == '>') 
		{
			if( sscanf(buf+1, "%s", next) != 1 ) {
				printf("bad format in %s or %s\n", buf, next);
				return 1;
			}
			else 
			{
				get_chr_name(next, chr);
			}
			printf(">%s\n", chr);
		}
		else {
			len = strlen(buf);
			if( buf[len-1] == '\n' ) printf("%s", buf);
			else printf("%s\n", buf);
		}
	}

	fclose(f);
	return EXIT_SUCCESS;
}
