include ../Flags.mk
DEST=../../bin

PROGS = dna seq_len gff_subset check_Ns common_reads fastq_head separate_fq edit_fq_head filter_protein_db_fasta check_splice_signal filter_gff comp_gffs check_aa dna2aa check_start_stop_seq check_start_stop_codons false_positive_orfs diff_list rm_redun_gff conv_scf_pos non_ref yeast_specific_common common_gene_list common_gene_set update_gff_blastx pull_c merge_gff gff2sim4 gff2codex sum_gff exonerate_gff2sim4 reverse_exon_order sum_all rm_endchar_fq trim_short_contigs  remove_lower_tri_and_diag conv_head extract_pos_info conv_scf_head readnames add_unique_only rm_redun_vcf make_table merge_bed vcf2gd_snp pull_fasta_scaf extract_scaf_codex diff_lines

all : $(PROGS)

dna:
	$(CC) $(CFLAGS) dna.c util.c seq.c nib.c charvec.c -o $(DEST)/dna

seq_len :
	$(CC) $(CFLAGS) seq_len.c util.c seq.c nib.c charvec.c -o $(DEST)/seq_len

gff_subset :
	$(CC) $(CFLAGS) gff_subset.c util_i.c util.c -o $(DEST)/gff_subset
	
check_Ns :
	$(CC) $(CFLAGS) check_Ns.c util.c seq.c nib.c charvec.c -o $(DEST)/check_Ns

common_reads :
	$(CC) $(CFLAGS) common_reads.c util.c -o $(DEST)/common_reads

fastq_head :
	$(CC) $(CFLAGS) fastq_head.c -o $(DEST)/fastq_head

separate_fq :
	$(CC) $(CFLAGS) separate_fq.c -o $(DEST)/separate_fq

edit_fq_head :
	$(CC) $(CFLAGS) edit_fq_head.c -o $(DEST)/edit_fq_head

filter_protein_db_fasta :
	$(CC) $(CFLAGS) filter_protein_db_fasta.c util.c -o $(DEST)/filter_protein_db_fasta

check_splice_signal :
	$(CC) $(CFLAGS) check_splice_signal.c tokens.c util_genes.c util_i.c util.c seq.c nib.c charvec.c -o $(DEST)/check_splice_signal


filter_gff :
	$(CC) $(CFLAGS) filter_gff.c tokens.c util_genes.c util_i.c util.c -o $(DEST)/filter_gff

comp_gffs :
	$(CC) $(CFLAGS) comp_gffs.c tokens.c util_genes.c util_i.c util.c -o $(DEST)/comp_gffs

check_aa:
	$(CC) $(CFLAGS) check_aa.c tokens.c util_genes.c util_i.c util.c -o $(DEST)/check_aa

dna2aa :
	$(CC) $(CFLAGS) dna2aa.c util.c -o $(DEST)/dna2aa

check_start_stop_seq :
	$(CC) $(CFLAGS) check_start_stop_seq.c util.c seq.c nib.c charvec.c -o $(DEST)/check_start_stop_seq

check_start_stop_codons :
	$(CC) $(CFLAGS) check_start_stop_codons.c util.c seq.c nib.c charvec.c -o $(DEST)/check_start_stop_codons

false_positive_orfs:
	$(CC) $(CFLAGS) false_positive_orfs.c util.c tokens.c util_i.c util_genes.c -o $(DEST)/false_positive_orfs

diff_list:
	$(CC) $(CFLAGS) diff_list.c util.c -o $(DEST)/diff_list

rm_redun_gff:
	$(CC) $(CFLAGS) rm_redun_gff.c tokens.c util_genes.c util_i.c util.c -o $(DEST)/rm_redun_gff

conv_scf_pos:
	$(CC) $(CFLAGS) conv_scf_pos.c util.c tokens.c util_i.c util_genes.c -o $(DEST)/conv_scf_pos

non_ref:
	$(CC) $(CFLAGS) non_ref.c util.c tokens.c util_i.c util_genes.c -o $(DEST)/non_ref

yeast_specific_common:
	$(CC) $(CFLAGS) yeast_specific_common.c util.c -o $(DEST)/yeast_specific_common

common_gene_list:
	$(CC) $(CFLAGS) common_gene_list.c util.c -o $(DEST)/common_gene_list

common_gene_set:
	$(CC) $(CFLAGS) common_gene_set.c util.c tokens.c util_i.c util_genes.c -o $(DEST)/common_gene_set

update_gff_blastx:
	$(CC) $(CFLAGS) update_gff_with_blastx_output.c util.c util_genes.c util_i.c tokens.c -o $(DEST)/update_gff_blastx	

pull_c :
	$(CC) $(CFLAGS) pull_c.c util.c seq.c nib.c charvec.c $(LFLAGS) -o $(DEST)/pull_c

merge_gff :
	$(CC) $(CFLAGS) merge_gff.c tokens.c util_genes.c util_i.c util.c -o $(DEST)/merge_gff

gff2sim4 :
	$(CC) $(CFLAGS) gff2sim4.c util.c -o $(DEST)/gff2sim4

gff2codex :
	$(CC) $(CFLAGS) gff2codex.c util.c -o $(DEST)/gff2codex

sum_gff :
	$(CC) $(CFLAGS) sum_gff.c util.c -o $(DEST)/sum_gff

exonerate_gff2sim4 :
	$(CC) $(CFLAGS) exonerate_gff2sim4.c util.c -o $(DEST)/exonerate_gff2sim4

reverse_exon_order :
	$(CC) $(CFLAGS) reverse_exon_order.c util.c -o $(DEST)/reverse_exon_order

sum_all :
	$(CC) $(CFLAGS) sum_all.c util.c -o $(DEST)/sum_all

rm_endchar_fq :
	$(CC) $(CFLAGS) rm_endchar_fq.c -o $(DEST)/rm_endchar_fq

trim_short_contigs : 
	$(CC) $(CFLAGS) trim_short_contigs.c -o $(DEST)/trim_short_contigs

remove_lower_tri_and_diag :
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c remove_lower_tri_and_diag.c -o $(DEST)/remove_lower_tri_and_diag

conv_head :
	$(CC) $(CFLAGS) conv_head.c -o $(DEST)/conv_head

extract_pos_info :
	$(CC) $(CFLAGS) extract_pos_info.c -o $(DEST)/extract_pos_info

conv_scf_head :
	$(CC) $(CFLAGS) conv_scf_head.c util.c -o $(DEST)/conv_scf_head

readnames :
	$(CC) $(CFLAGS) readnames.c -o $(DEST)/readnames

add_unique_only :
	$(CC) $(CFLAGS) util.c add_unique_only.c -o $(DEST)/add_unique_only

rm_redun_vcf :
	$(CC) $(CFLAGS) util.c rm_redun_vcf.c -o $(DEST)/rm_redun_vcf

make_table :
	$(CC) $(CFLAGS) make_table.c util.c -o $(DEST)/make_table

merge_bed :
	$(CC) $(CFLAGS) util.c util_i.c bed.c merge_bed.c -o $(DEST)/merge_bed

vcf2gd_snp :
	$(CC) $(CFLAGS) util.c util_i.c vcf2gd_snp.c -o $(DEST)/vcf2gd_snp

pull_fasta_scaf :
	$(CC) $(CFLAGS) pull_scaf_fasta.c -o $(DEST)/pull_fasta_scaf

extract_scaf_codex:
	$(CC) $(CFLAGS) extract_scaf_codex.c util.c -o $(DEST)/extract_scaf_codex

diff_lines:
	$(CC) $(CFLAGS) diff_lines.c util.c -o $(DEST)/diff_lines
