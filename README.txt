AGAPE Usage

Data Preparation

For a raw sequence in FASTQ,

1. Place the reference sequence file in FASTA and annotations in GFF in directory "reference" of the package directory (see file yeast.fasta and yeast.gff for file format)

2. Save the expressed sequencing tags (EST) and protein datasets for the reference in directory "cfg_files" of the package directory (see yeast_est.fasta and yeast_protein.fasta) and change the path of each file in configuration file "configs.cf" of the pacakge directory (REF_NAME, PROTEIN1, EST1, PROTEIN2, EST2, REF_DIR, REF_FASTA, and REF_GFF)

3. Creat a subdirectory in the package directory. 

4. Type "run_AGAPE_per_seq.sh Output_directory Output_prefix Fastq1 (or Fastq2)" to run all steps of the pipeline at once for assembly, annotation, and identification of novel genesi (note put the "full path" for the Fastq files). 

5. If you wish to run each step of the pipeline individually, do as follows.
	1) Assembly
	- type $(AGAPE_DIR)/agape_assembly.sh Output_directory Output_prefix AGAPE_directory Sequence_file1 (Sequence_file2 for paired-end)
	- output assembly file is saved in $(Output_directory)/$(Output_prefix).scf.fasta
	- Sequence_file1 and Sequence_file2 should include the full path of the file

	2) Annotation
	- type $(AGAPE_DIR)/agape_annot.sh Output_directory Output_prefix Assembly_FASTA_file AGAPE_directory
	- output annotation file is saved as
	  $(Output_directory)/comb_annot/$(Output_prefix).gff


	3) Identification of novel genes
	- type $(AGAPE_DIR)/agape_novel_genes.sh Output_directory Output_prefix AGAPE_directory Assembly_FASTA_file Annotation_GFF_file Fastq1 (or Fastq2 for paired end)
	- results file in GFF is $(Output_directory)/non_ref/$(Output_prefix).orfs.gff
