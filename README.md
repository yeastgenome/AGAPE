AGAPE
=====

AGAPE (Automated Genome Analysis PipelinE) for yeast pan-genome analysis is designed to automate the process of pan-genome analysis and encompasses assembly, annotation, and variation-calling steps.

-- INSTALLATION --

Prerequisite programs that need to get installed

(A) Assembly

1. ABYSS (http://www.bcgsc.ca/platform/bioinfo/software/abyss)
  - add "bin" directory of ABYSS to your $PATH
(e.g. if ABYSS bin directory is located in /home/gsong/AGAPE/programs/abyss-1.5.2/bin,
then edit line "export PATH=$PATH" to "PATH=/home/gsong/AGAPE/programs/abyss-1.5.2/bin:$PATH" in .bash_profile and reload it by running "source .bash_profile)
  - edit $ABYSS path in configs.cf file /home/AGAPE/programs/abyss-1.5.2/bin

2. SGA (https://github.com/jts/sga)
  - download SGA (e.g. typing 'wget https://github.com/jts/sga/archive/master.zipi')
  - SGA dependencies (google sparse library, bamtools library, and zlib)
  - ./autogen.sh
  - ./configure --with-sparsehash=/home/gsong/sparsehash --with-bamtools=/home/gsong/bamtools --prefix=/home/gsong/AGAPE/programs/SGA  (note bamtools should be later version than bamtools.2.3.0)
  - set PATH=/home/gsong/SGA/bin:$PATH
  - edit $SGA and $SGA_src paths in configs.cf

(B) Annotation

1. Reference genome sequence in FASTA and its annotations in GFF. If the reference genome and annotations are not available, place a sequence and annotations that are closest to your sequences
  - put the reference sequence to directory "reference" in the AGAPE main
    directory
  - specify the reference sequence file in $REF_FASTA and the reference
    annotation in $REF_GFF in the configs.cf file

2. Protein sequence database and EST sequence database for MAKER
  - store protein and EST database files in FASTA format in cfg_files
  - specify the two files names as $PROTEIN1 and $EST1 in configs.cf
  - default is yeast_protein.fasta and yeast_est.fasta in cfg_files
  - store another larger set of protein and EST database files than $PROTEIN1 and $EST1 and specify the larger files as $PROTEIN2 and $EST2 in configs.cf

3. NCBI-BLAST command line package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - BLAST depedency (boost library)
  - ./configure --with-boost=/home/gsong/boost --prefix=/home/gsong/AGAPE/programs/BLAST
  - type 'make'
  - type 'make install'
  - edit $BLAST path to '/home/gsong/programs/BLAST/bin' in configs.cf

3. ChainNet
  - Download executable program called 'axtChainNet.zip' in http://hgwdev.cse.ucsc.edu/~kent/exe/
  - put the ChainNet executable files to /home/gsong/AGAPE/programs/axtChainNet
  - check $axtChainNet path in configs.cf

4. Augustus (http://augustus.gobics.de)
  - set $AUGUSTUS path in configs.cf
  - set $AUGUSTUS_REF: the reference genome used for the AUGUSTUS annotation
    (see SPECIES.list file and pick the closest species of your genome
sequences)

5. MAKER
  1) maker_exe.ctl
  - edit the path of all executable programs that ar required for running MAKER

6. LASTZ (included in the pipeline)

(C) Variation calling

1. HUGESEQ

-- USAGE --
