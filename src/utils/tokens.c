#include "main.h"
#include "tokens.h"

void parse_last_column_gff(char *line, char *gname, char *seq_name, char *chr, int *b, int *e, float *pid)
{
	char chr_b[LEN_NAME];
	char chr_e[LEN_NAME];
	char pid_char[LEN_NAME];
	char temp1[LEN_NAME];
	char temp2[LEN_NAME];

	strcpy(chr_b, "");
	strcpy(chr_e, "");
	strcpy(pid_char, "");
	strcpy(temp1, "");
	strcpy(temp2, "");
	if( sscanf(line, "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%s", gname, seq_name, chr, chr_b, chr_e, temp1, temp2, pid_char) != 8 ) {
		fatalf("9th column format error: %s\n", line);
	}
	else {
		*b = atoi(chr_b);
		*e = atoi(chr_e);
		*pid = atof(pid_char);
	}
}

void parse_blastx_line(char *line, char *info, int *b, int *e, int type)
{
	char evalue[LEN_NAME], name[LEN_NAME];
	char gene_name[LEN_NAME];
	char pos[LEN_NAME];
	char cur_pos[LEN_NAME];
	char chr[LEN_NAME];
	char seq_name[LEN_NAME];
	char assem_type[LEN_NAME];
	char char_b[LEN_NAME], char_e[LEN_NAME];
	int len = 0;
	int i = 0, j = 0;
	int cur_b = 0, cur_e = 0;
	int num_exons = 0;
	
	strcpy(seq_name, "");
	strcpy(assem_type, "");
	strcpy(gene_name, "");
	strcpy(name, "");
	strcpy(pos, "");
	strcpy(cur_pos, "");
	strcpy(evalue, "");
	strcpy(chr, "");
	strcpy(char_b, "");
	strcpy(char_e, "");

	if( type == SGD ) {
		if( sscanf(line, "%s %*s %s %s %*s %*s %s %*s %s %*s", name, evalue, gene_name, chr, pos) != 5 ) {
			fatalf("wrong blast line in %s", line);
		}
		len = strlen(pos);
		i = 0;
		while(i < len) {
			if( pos[i] != ',' ) {
				cur_pos[j] = pos[i];
				j++;
			}
			else {
				cur_pos[j] = '\0';
				split_b_and_e(cur_pos, &cur_b, &cur_e);
				if( num_exons == 0 ) *b = cur_b;
				num_exons++;
				j = 0;
			}
			i++;
		}	
		*e = cur_e;
		sprintf(info, "SacCer,chr%s,%d,%d,%s,%s", chr, *b, *e, gene_name, evalue);
	}
	else if( type == ENSEMBL ) {
		if( sscanf(line, "%s %*s %s %*s %s %*s", gene_name, evalue, name) != 3 ) {
			fatalf("wrong blast line in %s", line);
		}
		else {
			if( sscanf(name, "%[^:]:%[^:]:%[^:]:%[^:]:%[^:]:%s", assem_type, seq_name, chr, char_b, char_e, pos) != 6 ) {
				fatalf("wrong blast line with ENSEMBL DB in %s", name);
			}
			else {
				strcpy(pos, char_b);
				strcat(pos, "-");
				strcat(pos, char_e);
				split_b_and_e(pos, b, e);
				
				if( strstr(assem_type, "chrom") != NULL ) {
					if( strstr(seq_name, "SacCer") != NULL ) {
						sprintf(info, "SacCer,chr%s,%d,%d,%s,%s", chr, *b, *e, gene_name, evalue);
					}
					else {
						sprintf(info, "%s,chr%s,%d,%d,%s,%s", seq_name, chr, *b, *e, gene_name, evalue);
					}
				}
				else if( strstr(assem_type, "contig") != NULL ) {
						sprintf(info, "%s,%s,%d,%d,%s,%s", seq_name, chr, *b, *e, gene_name, evalue);
				}
				else {
						sprintf(info, "%s,%s,%d,%d,%s,%s", seq_name, chr, *b, *e, gene_name, evalue);
				}
			}
		}

	}
	else {
		fatalf("unsupportive type: %d\n", type); 
	}
}

void parse_line(char *line)
{
	int i = 0;
	char index[LEN_NAME], name[LEN_NAME];
	char gene_name[LEN_NAME];
	char str = '\0';
	int len = 0;
	
	len = strlen(line);
	gene_name[0] = '\0';
	index[0] = '\0';
	name[0] = '\0';

	str = line[0];
	while( (i < len) && ( line[i] != '\0' ) && (line[i] != '\n') ) {
		i = concat_tokens_on_bracket(line, i, index, name);
		if( strcmp(index, "gene") == 0 ) {
			strcpy(gene_name, name);
		}

		if( strcmp(index, "location") == 0 ) {
			print_gene(str, gene_name, name);
		}
	}
}

void print_gene(char str, char *gname, char *line)
{
	int *b, *e;

	b = (int *) ckalloc(sizeof(int));
	e = (int *) ckalloc(sizeof(int));

	get_start_and_end(line, b, e);

	printf("%c %d %d %s\n", str, *b, *e, gname);
	print_exons(line);

	free(b);
	free(e);
}

void rm_temp_num(char *token, char *name)
{
	int i = 0;
	int len = 0;

	len = strlen(token);
	i = len-1;
	strcpy(name, token);
	while( (i >= 0) && (token[i] != '-') ) i--;

	if( token[i] == '-' ) {
		name[i] = '\0';
	}
	else {
		fatalf("field %s does not contain seperator -\n", token);
	}
}

void split_b_and_e(char *token, int *b, int *e)
{
	int i = 0;
	int len = 0;
	char *num;
	int temp = 0;

	len = strlen(token);
	num = (char *) ckalloc(len * sizeof(char));
	while( ((*token) != '\0') && ((*token) != '\n') && ((*token) != '-')) {
		num[i] = *token;
		i++;
		token++;
	}
	if( (*token) == '-' ) {
		token++;
	}
	else {
		fatalf("field %s does not contain seperator -\n", token);
	}
	num[i] = '\0';	
	*b = atoi(num);
	*e = atoi(token);	
		
	if( *e < *b ) {
		temp = *e;
		*e = *b;
		*b = temp;
	}
	sprintf(token, "%d-%d", *b, *e);
	free(num);
}

void get_start_and_end(char *line, int *b, int *e)
{
	int len = 0;
	int i = 0, j = 0;
	char num[LEN_NAME];

	*b = -1;
	*e = -1;
	num[0] = '\0';
	len = strlen(line);
  if( line[i] != '(' ) {
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != '(')) i++;
  }

  if( line[i] == '(' ) {
    i++;
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != '.') )
    {
      num[j] = line[i];
      j++;
      if( j >= LEN_NAME ) {
        fatalf("Over the max length for a gene name ( 100 characters ) in %s\n", line);
      }
      i++;
    }

		num[j] = '\0';
		*b = atoi(num);

    if( line[i] == '.' ) {
			i = i+2;
			j = 0;
    	while((line[i] != '\0') && (line[i] != '\n') && (line[i] != ')') )
    	{
     		num[j] = line[i];
      	j++;
				if( line[i] == '.' ) j = 0;
				i++;
			}
			num[j] = '\0';
		}
		*e = atoi(num);
	}
}

void print_exons(char *line)
{
  int len = 0;
  int i = 0, j = 0;
  char num[LEN_NAME];
	int num1 = 0, num2 = 0;

  num[0] = '\0';
  len = strlen(line);
  if( line[i] != '(' ) {
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != '(')) i++;
  }

  if( line[i] == '(' ) {
    i++;
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != ')') )
		{
    	while((line[i] != '\0') && (line[i] != '\n') && (line[i] != ',') && (line[i] != ')'))    
			{
     		num[j] = line[i];
      	j++;      
				if( j >= LEN_NAME ) {
        	fatalf("Over the max length for a gene name ( 100 characters ) in %s\n", line);
      	}
				if( line[i] == '.' ) {
					num[j] = '\0';
					num1 = atoi(num);
					j = 0;
					i++;
				}
      	i++;
			}
			if( (line[i] == ',') || (line[i] == ')') ) {
				num[j] = '\0';		
				num2 = atoi(num);
				printf("%d %d\n", num1, num2);
				j = 0;
				i++;
			}
    }
  }
}

int concat_tokens_on_bracket(char *line, int loc, char *index, char *name)
{
  int len = 0, i = 0, j = 0, k = 0;

  len = strlen(line);
  i = loc;
  if( line[i] != '[' ) {
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != '[')) i++;
  }

  if( line[i] == '[' ) {
    i++;
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] != '=') )
    {
      index[j] = line[i];
      j++;
      if( j >= LEN_NAME ) {
        fatalf("Over the max length for a gene name ( 100 characters ) in %s\n", line);
      }
      i++;
    }

    if( line[i] == '=' ) {
      index[j] = '\0';
      i++;

      while((line[i] != '\0') && (line[i] != '\n') && (line[i] != ']') )
      {
        name[k] = line[i];
        k++;
        if( k >= LEN_NAME ) {
          fatalf("Over the max length for a gene name ( 100 characters ) in %s\n", line);
        }
        i++;
      }
			name[k] = '\0';
    }
    else {
      name[k] = '\0';
    }
  }
  else {
    index[j] = '\0';
    name[k] = '\0';
  }

  return(i);
}

double AtoF(char *s)
{
	double val = (double) 0.0, power = (double) 1.0;
	int i = 0, sign = -1, exp = 0;
    
	for(i = 0; isspace(s[i]); i++) //skip white space
	;
	sign = (s[i] == '-') ? -1 : 1;
	if(s[i] == '+' || s[i] == '-')
		i++;
	for(val = 0.0; isdigit(s[i]); i++)
		val = 10.0 * val +(s[i] -'0');
	if(s[i] == '.')
		i++;
	for(power = 1.0; isdigit(s[i]); i++)
	{
		val = 10.0 * val + (s[i] -'0');
		power *= 10.0;
	}
	val = sign * val / power;
	if(s[i] == 'e' || s[i] == 'E')
	{
		i++;
		sign = (s[i] == '-') ? -1 : 1;
		if(s[i] == '+' || s[i] == '-')
		i++;
		for(exp = 0; isdigit(s[i]); i++)
		exp = 10*exp + (s[i] - '0');
		if(sign == 1)
		{
			for(; exp > 0; exp--)
				val *=10;
		}
		else
		{
			for(; exp > 0; exp--)
				val /= 10;
		}
	}
	return val;
}
