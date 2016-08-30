This file is the user manual for ISQuest_1.4.1...

1. Building ISQuest Executable
    a. Download and unzip the code in a directory.
    b. Change to the directory.
    c. Execute command: make

2. Installation

The ISQuest tool requires some local setup of supporting tools and databases.   

2a. Installation on Windows
	a. ISQuest executable file generated in step 1. 
	b. Download and install Microsoft Visual C++ 2010 Redistributable Package (x64)
		http://www.microsoft.com/en-us/download/confirmation.aspx?id=14632
	c. Install standalone BLAST. See manual 
		http://www.ncbi.nlm.nih.gov/books/NBK52637/
	d. Download "nt" database and setup as instructed in step(c).
	e. Edit PATH environment variable to include the path to isq executable file.

2b. Installation on Linux
	a. ISQuest executable file generated in step 1. 
	b. Install standalone BLAST. See manual http://www.ncbi.nlm.nih.gov/books/NBK52637/
	c. Download "nt" database and setup as instructed in step(c).
	d. Edit PATH environment variable to include the path to isq executable file.

2b. Installation with MPI
	a. mpic++ should be in PATH.
	b. Run make as usual. A new executable named mpiblast will be built.
	c. Path to mpiblast must be available in the config file as MPIBLAST. This it to make sure this can be invoked from SGE.
	d. Set MPIHOSTCOUNT can be set. Default is 4. 
	e. Run isq as "qsub -pe <env> 1 isq ...."


3. User Instructions

To use the software you will need:

	a. Parameter file with the following fields:
	    GENBANKPATH=<<SOME PATH>>
		BLASTPATH=<SOME PATH>\blast-2.2.29\bin\
		BLASTDB=<<SOME PATH>>\blastdb\
		LOOPCOUNTER=2
		USEPREVBLAST=1
		MPIBLAST=/scratch/abisw001/ISQuest/isquest-code/
		MPIHOSTCOUNT=8
	b. A set of contigs in fasta format.
	c. The set of reads used to assemble the contigs in fasta format.
	d. Prefix of the name of the output files.

The command to execute the file is as follows:
	
	>isq <paramter file> <output directory> <contig file> <read file> strain_name
	
4. Important Points
	a. The software during the first run will download a lot of GenBank files in GENBANKPATH from NCBI. 
	   This will make this run very slow. Future runs will be faster. 
	b. The software will do a BLAST search during the first run. This will make the first run slow. 
       Future runs can use the same BLAST results. If you want to generate new BLAST search use 
 	   USEPREVBLAST=0
	c. The GBFP parser that I am using right now cannot handle certain GenBank file types. 
	   There is no support for this anymore and I have to go in and fix the bugs.
	   This will be done in the next release. 


5. Output Files
	a. Bunch of *.sam alignment files with corresponding fasta files with same name containing the reference sequence 
		i. You can import them into Geneious
		ii. Naming Convention <strain_name>_<loopID>_transposaseHit_<level 1 ID>_<level 2 ID>.*.sam
	b. IS_Contig_Map.txt provides the contig ID and the coordinates of the original transposase hits. 
