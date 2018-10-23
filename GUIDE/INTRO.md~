#Annocript 0.2 - Introduction

Annocript is a pipeline that can be used to manage and give a consistence in terms of annotation to de-novo generated transcriptomes. It executes several blast analysis steps and builds annotations using information taken from UniProt, NCBI Conserved Domain Database and Nucleotide divisions (CDD), Gene Ontology, and the Enzyme Commission. *Annocript* also perform analyses based on Portrait software to measure the noncoding potential of RNA sequences. *Annocript* works with both nucleotides and peptides. If you use proteins different blast programs will be ran while Portrait and dna2pep will not, since they work only with nucleotides.

Notice: *Annocript* has been tested only on Linux Ubuntu 64 bit Systems. We are not sure of what will happen on different operating systems and machines!!

##Methods

The pipeline allows the creation of a comprehensive user-friendly table containing all the annotations produced for each transcript. The user can choose to annotate her/his transcriptome against selected organisms or the complete database. 

 The following are the analyses performed:
- The proteins most similar to the transcripts are given by the blastx (blastp if you use peptides) analyses against the UniProt databases SwissProt and TrEMBL (or UniRef). Default databases are downloaded at:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note

- Blastn (tblastn) against a concatenation of the SILVA database (small and large subunits ribosomal RNAs) and the Rfam database allows to check for ribosomal and other short noncoding RNAs. 
Default databases are downloaded at:
  - ftp://ftp.arb-silva.de/current/Exports/SILVA_119_LSUParc_tax_silva.fasta.gz
  - ftp://ftp.arb-silva.de/current/Exports/SILVA_119_SSUParc_tax_silva.fasta.gz
  - ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz

- Rpstblastn (rpsblast) returns information about the Conserved Domains Database within each transcript. 
Default database is downloaded at:
  - ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz

- Mapping of GO functional classification is shown using the best matches between SwissProt and TrEMBL. If UniRef is used, the GO terms are always taken associated to its result. GO terms can be also associated to Pfam Domains. To choose the association type (protein,domains or both) you can use the configuration file.
The mapping of GO terms with UniProt identifiers is obtained using the file at:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

- Mapping of Enzyme Commission IDs and Pathways descriptions are always given associated only to the SwissProt id, if present. The mapping of these identifiers with UniProt ones is obtained using the files at:
  - ftp://ftp.expasy.org/databases/enzyme/enzyme.dat
  - http://www.uniprot.org/docs/pathway.txt


- Portrait measures the probability that a sequence is coding or non-coding and its score, together with a final heuristic, based on the integration of all the results, makes Annocript capable to also identify bona-fide noncoding transcripts.

- An HTML page is generated containing summary plots and descriptive statistics of a given transcriptome. GFF3 outputs allow the integration of the annotations in genome browsers and easy file exchange trough current standard formats and a GFF database permits to rapidly access results using bioinformatics API such as BIOPERL.

Let's use it!

When you decide to use *Annocript* you can choose to run one or all of the following 4 steps.

- Database creation: it creates a large database with informations about proteins from UniProt (creation of the source annotation database is mandatory);
- Execution of programs: it runs the blast programs, Portrait to measure the noncoding potential and dna2pep to search for the longest ORFs;
- Output creation: it creates a large file with the programs results and create a GFF database and GFF3 files;
 Creation of statistics: this creates a web page with statistics. 

You should create the database only once and then *Annocript* will use it. (You don't need to repeat the creation a second time unless UniProt has produced another version of UniProt knowledge base and you want to update it). You should then run the programs to obtain the analysis output. To customize the execution you have to use the configuration file.

##Configuration Files

*Annocript* contains two different configuration files. One is for users and the other is used by *Annocript* to choose names and other stuff that should not be changed by the generic *Annocript* user. However, experienced power-users can see how to use the latter ADVANCED_USAGE section of the Guide. By modifying the parameters in the user configuration file you can run only some programs or only have a database for doing analysis using your scripts, you can set the parameters and the cutoff for the analysis and so on.
The user configuration file is divided in 3 main parts. They are the following:

- the first one is composed by simple flags to indicate which programs to run and which not;
- the second part is composed of parameters the inputs to be used in the programs used by Annocript;
- the last part is with paths to the programs and to the online databases.

Examples of whole transcriptome annotations made by Annocript can be downloaded from http://bit.ly/15vnALW. 
---------------------------------


If you get some error during the installation or the running of Annocript please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
