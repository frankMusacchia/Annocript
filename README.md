#*Annocript*
==========

**Annocript** is a pipeline for the annotation of *de-novo* generated transcriptomes. It executes several BLAST
 analysis steps with UniProt, NCBI Conserved Domain Database and Nucleotide divisions adding also annotations 
from Gene Ontology, the Enzyme Commission and Pathways. *Annocript* also uses PORTRAIT software to measure non-coding 
potential.

##Mission
*Annocript* is a highly configurable tool to annotate and analyze transcriptomes. Since it is not coming with a
 double-click installation you should dedicate a bit of time the first time you use it. But don't run away! 
This is a very simple guide, it is big only because everything is explained step-by-step but we promise that,
 after you installed and ran once, you will do it very rapidly the following times and you will never stop using it!
We strongly suggest to keep open this guide while you are installing and running *Annocript* for the first time.

##System Requirements

To work, **Annocript** needs:
- MySQL database server (tested with ver5.1) 
- You also need a MySQL client (tested with ver5.5, with a working user account with user-id and password)
- PERL (tested with ver5.10) 
- BIOPERL (tested with ver1.6)
- PYTHON (tested with ver2.6)
- R (tested with ver3.0)
- BLAST+ (tested with ver2.25)
- PORTRAIT (tested with ver1.1)
- DNA2PEP (tested with ver1.1)


##Installation and tutorial

To start using it you must follow instructions in the **GUIDE** folder.
Below is a simple introduction of what you will find there and how to find solutions to possible problems.
- **Annocript** has some specific purposes. To figure out what they are, please read the **INTRO.txt** file;
- Since *Annocript* is a pipeline, it uses some software that you have to download and install. Please follow 
the instructions on how to do it in the **INSTALL.txt** file. In there you will also learn how to install *Annocript*;
- Every software has is thorns... To understand how to run *Annocript* the first time please read the TUTORIAL.txt file;
- When you will be an advanced user you will want to use *Annocript* in a faster and smarter way. ADVANCED_USAGE.txt 
explains how to begin to do this. Here, an HINT and TIPS section contains a lot of useful information on how to
 perform advanced operations;
- **FAQ.txt** contains Frequently Asked Questions (refer to the FAQ whenever you have errors that block *Annocript*);
- **OUTPUT.TXT** contains a list of all the output that *Annocript* gives, the organization of folders and 
of the configuration files;
- **HISTORY.txt** is the history of changes to *Annocript* from it's creation to the current release.


Please enjoy this product and send us questions and suggestions!


##*Annocript* times and space

Any user will wonder how much time and space *Annocript* needs. We can say approximate times and sure space dimensions.
*Annocript* execution time depends strongly from what modules the user is executing and of course from the speed and memory
of your machine. Here we only show the most important processes that are always taking the biggest part of the overall 
computational time.


###Database Creation

The database creation needs to be done only once but it depends on what proteins db you are creating (Uniref or Uniprot_kb).
If you decide to use Uniprot_kb, *Annocript* will fill an annotation table with both the BLAST results from Swiss-Prot and Trembl
while if you are using Uniref, both Swiss-Prot and Uniref result will be shown.
The database comprise also information about domains, GO, enzymes and pathways. Furthermore, databases to use for BLAST will 
be downloaded also (see the **INTRO** file for more information about this step).

Consider as a suggestion that, if your machine is not very recent and powerful you can consider to install Uniref database.
Otherwise Uniprot_kb is a good choice.

####Times

We took times for the creation of the database from 3 different machines with different cores and RAM memory: 
- Machine with 24 core and 96 GBRAM 
	- Uniprot_kb_2014_03 - 4 hours
	- Uniref_2014_04 - 3 hours
- Machine with 24 core and 24 GBRAM 
	- Uniref 2013_01 - 9 hours 
	- Uniprot_kb 2014_03 - 21 hours 
- Machine with 8 core and 8 GBRAM
	- Uniref 2014_05 - 38 hours

The speed of the overall process depends greatly from the time required for the download of the huge databases. Here the times were 
taken by using approximately 1MB/sec in download but depending from the speed of you connection the download process may 
take few hours more.


####Space

The space needed to build the MySQL database is approximately: 
- UNIREF 2014_04: 52.61GB (50.1GB for files; 2.4 GB MySQL db) 
- UNIPROTKB 2014_04: 92.6GB (85.3GB for files; 7.3 GB MySQL db)

Newer protein databases will slightly increase the necessary space according to their growth in term of proteins numbers.
You can save some space by choosing, during the configuration, to remove the downloaded .zip files.
- zip files in Uniref: 8.5GB
- zip files in Uniprotkb: 14.5GB

###Programs Execution

Programs execution strongly depends from the transcriptome and the database we are using. The number of sequences but mostly 
their mean length.
We took the times of the BLASTx execution (the most expensive process) against the Uniref90 database with our default parameters 
(which you will find in the CONFIG_USER.txt) on a machine with 24 threads (2 Intel(R) Xeon(R) X5660 at 2.80GHz and 96GB RAM). 
This machine permitted the multi-threaded run of blast and rpsblast. We always used 20 core. The followings are the times when the
 given mean length of transcriptome.

####Times

It needed about: 
- 6 hrs for 11351 transcripts of mean length 951 nucleotides; 
- 13 hrs for 30346 transcrips of mean length 1245 nucleotides;
- 20 hrs for 17776 transcripts of mean length 1478 nucleotides;
- 23 hrs for 30339 transcripts of mean length 1427 nucleotides;
- 46 hrs for 98174 transcripts of mean length 965 nucleotides;
- 72 hrs for 64388 transcripts of mean length 2125 nucleotides.

Depending by the powerful of you machine, the BLASTx and the rpsBLAST execution may take few hours more. 


Find me at: francesco [dot] musacchia [at] szn [dot]it