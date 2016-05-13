#The *Annocript* 0.2 - History

*Annocript* 0.1 (2011) - Completely realized by Swaraj Basu under the supervision of Remo Sanges. It was based on the manual download of the databases and on the use of a configuration file that was created at running time at the first execution of the software. Then *Annocript* creates a little database with information coming from the output of BLAST software.

*Annocript* 0.2. In version 0.2 the software was completely rewritten by Francesco Musacchia under the supervision of Remo Sanges. It creates a database with informations about proteins by downloading:
- UniProt Knowledge Base and merging the Swiss-Prot and TrEMBL datasets;
- or UniRef90 and adding the results of Swiss-Prot too.
The code for the execution of the BLAST software is copied from Annocript 0.1

All the *Annocript* has been divided in four main modules:
- the DB creation, that is needed for the creation of a large database with informations taken from UniProt, UniRef... etc, about proteins
- the programs execution, that is a script that permits the executions of the various program that we want to
use (BLASTx, BLASTn, Portrait, ..etc). Its purpose is to execute programs and store outputs.
- the output creation with GFF3 use, this is the part of *Annocript* where the programs output is extracted and 
inserted in a final output table. This output can also be converted in GFF3 files and stored with BioPerl;
- the creation of statistics is a script that extract some files with statistics from the output table and creates a web page with plots.


*Annocript* 0.2.10: it uses 3 folders with 3 different execution based on 3 differents config files. Here there's the problem that the output files have different name with the same files that areused from another execution as input.	Moreover the use of 3 file of configuration is too much heavy for usability of the software.

*Annocript* 0.2.11: it uses only 1 configuration file and then the first script is the executor that calls all other scripts and manages their config hashes. The programs have to save all the file in a folder data.		For each execution now the software asks for the folder of user session. This could be existent or not. It can be created and overwrited.	In this version I've also done some changes to the software to keep it more readable	and fast in example changing the database names. Here I also managed the fasta file to insert into the Portrait software using some scripts wrote by Giuseppe.

*Annocript* 0.2.12: I decided to make some changes to the database: the mapping table will have the id for each enzymeId since I was making two queries each time and the computational time was very high. This is exactly in the step where I build the final table in GFF3_AND_OUTPUT and the system does several queries to find what is the best between the trembl and sp result. The problem is how to reduce this computational time. So I started a new version.

*Annocript* 0.2.13: this version is characterized from the different production of the final table. It was all done by accessing the GFF database and it was very expensive. The time to only create the table was very high.
		So we decided to use other mechanisms and finally I designed an hash from where data are taken instead that make queries. This has been definitively more efficient. All the script for GFF producing and final table occurs 3-4 hours.	In this version has been established that there will be one only folder for all the data that is computed by the EXECUTE_PROGRAMS script when the blasts have to access the databases (created by makeblastdb and the rFam and CDD). This folder is located in the main folder of *Annocript*. (Fri, 07 december 2012)...*Annocript* is becoming a monster...

*Annocript* 0.2.14: this is a version that we begun to use in several projects. It was almost complete in all the initially designed functionalities and it contains possibility to compute some statistics on the output files and creates also some plots with R.	All these plots are put in a HTML file with also the statistics.

*Annocript* 0.2.15: in this version two additions have been done. The first is that you can have different files of output for the BLAST programs. In this way one can choose to divide the file of transcripts and run them in parts on different machines or processors. Another main addition is that now *Annocript* can be ran in background. We used a perl script that executes the configuration and another script that runs all the	pipeline in background.

*Annocript* 0.2.16: in this version I used only one configuration hash to memorize variables. Then with Remo we seen the code together and some adjustement were done.    

*Annocript* 0.2.17: this is a version wrote for the publication. A lot of IF has been removed to permit the continue 
 execution of all the pipeline without too questions to the user. Some other improvements are: to check if a database exists we used  the query 'show database like'; used two log file to know what elements of the GOMapping are not found in the UNIPROT database; the module gff3_utilities has been removed and the functions to have the ORF file have been put in the gff3_manager module; built the functions dlAndExtract to download and extract a compressed file; for the uncompression of files has been now used a native function of PERL that finds the compression type (subroutine extract_file); a new function is also extract name to extract the name of the file without extension (extractName);the two functions similar to create a NonCoding info file and to know what is the number of non coding sequences have been joined in one only;the query coverage and hit coverage have been added to the final table  for both swiss-prot and trembl results;

*Annocript* 0.2.18: in this version we changed the way the mappingTable is created. I was using an INSERT with the option to check if the element is already present it will be updated was too expensive in computational terms so we decided to use an hash built during the reading of the files. This results in a faster process since there is only one access to the database. Here, thanks to Giuseppe's feedback and suggestions, we added the possibility to use the UniRef database. So a new function for the parsing has been added (createUnirefDB).Now *Annocript* searches for both an UniRef and an UNIPROT_KB database and asks what to use. An error has been corrected in the output_and_stats module since it was counting the GO occurrences multiple times for the same sequence annotation.
A new field has been added to the final table saying if the sequence is coding or non-coding for *Annocript*. This is done by changing the printAllInfo subroutine in the module and by checking the constrains for each of the blastx, rps, blastn outputs and	the longest ORF lenght and the non coding probability given by portrait.	This resulted also in the possibility to change completely the function createNCInfoFile that was doing all these checks and now can simply look for the corresponding field to create the file.

*Annocript* 0.2.19: in this version I read the code with my supervisor Remo and we corrected some stuffs in 		output_and_stats.pm: the count of annotated sequences is done during the building the array of all results; the non-codingness of a sequence is computed during the construction of the final output table	since we need to put a new field with this flag; the function execute_stats has been modified since it was using equal pieces of code:  the subroutines extractGOstats, printOutGOStats and extractCDDstats have been created to reduce the amount of code repeated;

*Annocript* 0.2.20: in this version we put also the first part of the print in the LOG file because it was not there and can be useful for the user: functions annoPrint and annoDie have been introduced in the *Annocript*.pl
file for this. There was an error in the count of the enrichments. It has been solved. 	In this version also have been extracted some common use functions from the modules and inserted in a different module called utilities.pm.

*Annocript* 0.2.21: Here *Annocript* has been corrected after the alpha-test executed by Marco.

*Annocript* 0.2.22: This is a version with some other corrections. Now *Annocript* is ready for the beta-test. *Annocript* has also been presented to BITS2014 and has been awarded in the best poster competition.

*Annocript* 0.2.23: After the beta test some corrections had been done. In particular we had the problem with the download using LWP and with Head function to simply access a link and see if it is good. The problem was of the institute so I decided to insert an alert. Some corrections have been done in the user guide and little problem came out are all solved. Another problem we had with the use of the Extract file function that is having a bad behaviour. The uncompression on slower machines is often blocked. I decided to write my own function to uncompress the files. 
The problems in this part are mostly coming from the usage of *Annocript* on different machines. 

*Annocript* 0.2.23: We finished to correct the software and on 1th April 2014 I gave a presentation with all the issues had during the beta-test and Remo and Swaraj gave suggestion to implement a table where one can write the table if it has been well written or Nothing otherwise. If  it is there, the table will  not be built again in a different *Annocript* execution.
The second suggestion is to use the CUT command to cut the mapping file from Uniprot to catch only the GO ids. This should reduce the computational time

*Annocript* 0.2.24: The second *Annocript* day has been done but the process stopped at the loading of the database in less powerful systems. It was an error in system execution with a 'mysqlimport' not taking in account that the string interpolation was reading '\n' as newline and '\t' as a tab. In this new version different things are changed. 
	-the table upload in db is did first with DBI, thus a query to DB, and if this does not work (you can have an error: DBI error: The used command is not allowed with this MySQL version) a system call is doing the same job;
	-system calls now are did with a set of controls shown on the PERL reference website to be more
specific;
	-tables are not more in the *Annocript* folder while they are being generated. Instead they will be in the data/YEAR_MONTH folder in an apposite TABLE folder;
	-the config hash before was saved in the *Annocript* folder. Now it will stay in the SESSION/DATA folder and it is passed as argument to *Annocript*_executor.pl;
	-in *Annocript* folder now compare only few files. Now folders contain the useful files and script:
		- CONFIGURATION: contains everything regarding the configuration
		- USEFUL: contains the scripts utilities.pm, utils.pm and rRNA zip file

*Annocript* 0.2.25: *Annocript* till this version was using a lot of memory to create the mapping table. In this new version we decided to insert a control and if there is enough memory then *Annocript* will use a big hash to store information while indeed it will use queries to access the database.
In this version also has been completely removed the use of the file with Counts since a single file is not useful.

*Annocript* 0.2.26_DELIVERED: *Annocript* took is final shape in this version. New adjustments where did and tests were performed. It was approximately september 2014 when I decided to call this the first release version. We sent the paper in october using this last version

*Annocript* 0.2.28:  Before to know if the paper was published or not I did some little adjustments to the overall pipeline. In particular with the versione 0.2.7: 
- all the downloads are did as first step, in this way we are not continuosly dependent from the connection, once that the databases are downloaded *Annocript* has only to build its own db; 
- we added columns to the table to give also the strand corresponding to the BLAST alignment; 
- a parameter is added to decide how many elements should be given in the plots, before they were 20;

*Annocript* 0.2.28b: 
- Now the table of output in Annocript contains a column for each of the levels of pathways. 
- Added plots with statistics about lengths of the transcripts and coverage of the results obtained from blast; 
- removed every warning coming from the function delete_file in utilities.pm when files are not found 
- printed dna2pep command.
- I modified the code for the HTML page. Now we can construct more than one page and there is an indipendent module to write them. 
- Also the appearence of the HTML result is changed: since they are now three pages, I added a section with links for the navigation. In this version for the first time appears the logo of Annocript in the HTML results. 
- I added more specificity with the use of Uniref. Now Annocript can install only one Uniref database type per each version downloaded (50,90,100) and the uniref db type is searched in the folder and not in the configuration file. After an update the configuration file will change and will take the default (50).
- Added suggestion after each of the "strange" questions of Annocript. Like "should I remove the zip files?". Even if they are considered simple questions, they could be confusing for users because they do not know what to do.
- Now the folder with the slices of the fasta can be removed if the option "keepTempFiles" is set to NO (default behaviour).If the user sets YES, they will be kept.

*Annocript* 0.2.29
At the end of november 2014 we received the comments from the reviewers of Bioinformatics. They asked for two main modifies at the software. One was to show GO terms assigned to domains instead that proteins. I decided to add for doing both the annotation. Thus, the user can choose to assign GO terms to both proteins and domains or to one of them exclusively. Another hint received by the reviewers was to be possible to choose the organisms to BLAST against. Thus, I modified the code so that the user can give a file in input where he/she can write the names of the organisms and *Annocript* will filter the complete TrEMBL and SwissProt databases to obtain a reduced BLAST index.

*Annocript* 0.2.30

- Fixed a bug that did not import the database schema when imported on different host (hint by Miquel on the forum)
- Removed the customized database of ribosomal and adopted a new policy for building the database of noncoding. The database is downloaded from Silva and joined Rfam.
- Control of the FASTA file. If you do not meet these conditions, you are notified and Annocript stops: the name of the contigs always immediately after the>; Max 80 characters per line; only the presence of characters IUPAC.
- I added the possibility to do the BLASTx against only SwissProt or TrEMBL/UniRef. Instead the user can choose to get only the result from one of the two by shifting the corresponding blast ouput file in another folder.
- Removed the colorBrewer that caused trouble in the auto-installation. Inserted a manual staining.
- Fixed an error in the management of Portrait's output.
- Fixed the error that if there is not the result of dna2pep in the output table PERL exits the program.
- Added a check for the actual construction of the indexed database after makeblastdb with a function that checks the existence of all files usually created by the program