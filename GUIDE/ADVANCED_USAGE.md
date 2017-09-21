#Annocript 0.2 - Advanced usage

If you have already experienced several Annocript successful runs you may want to perform some more specific tasks or just using some sections of it (e.g. create only a database, running only some programs 
or only to obtain results). You can do it buy using the user configuration file in a very simple way. 
Please read the following tutorials to do it! 


##Build a new Protein database 


The first thing that you need for you analysis is to install a proteins database with UniRef or UniProt_KB (knowledge base) 
To build a database open the configuration file in the Annocript folder and set parameters as follows:

                ---------------- TERMINAL ------------------ 
                francesco@compaq2:~/ann_works$ nano config_user.txt 
                ---------------- TERMINAL ------------------ 

Set the name of your fasta sequences file (or a fake one) and choose to execute only the database creation: 

                ---------------- TERMINAL ------------------ 
                fastaSeqs = sequences.fasta 
                
                doDbCreation = YES 
                doExecutePrograms = NO 
                doBuildOutput = NO 
                extractStatistics = NO 
                ---------------- TERMINAL ------------------ 


Put the file with the fasta sequences in your working folder and set the name of the fasta file with sequences: 
fastaSeqs = yourFastaSequences.fasta 
(Since Annocript has been created to do analysis you must set  the fastaSeqs parameters. If you don't want to set 
because you need only a database, please leave the default name. Annocript is smart: if you don't choose to do analysis, 
it will not bother you saying that you don't have transcripts file.) 


On command line go to your chosen folder and run Annocript (the config_user.txt file must stay in your working directory) 

                ---------------- TERMINAL ------------------ 
                francesco@compaq2:~/ann_works$../Annocript/annocript.pl config_user.txt 
                
        
                1. New session - 2. Use Previous - 3. Exit 
                Your choose:  1 

                Choose a name for your session: new_database 
                /home/francesco/ann_works/jobs/new_database created! 

                Copying the given config file in /home/francesco/ann_works/jobs/new_database 
                ##################################################### 
                # CHECKING THE INTEGRITY OF THE CONFIGURATION FILES # 
                ##################################################### 
                They are good. Annocript can proceed! 
                ################################### 
                # RUNNING ANNOCRIPT CONFIGURATION # 
                ################################### 
                DB file will be downloaded from UNIPROT website and DB built. 
                The final table will be built 
                A web page and file with statistics will be created 
                Do you agree with it?(y or n) y 
                Ok!.. let's go! 
                ---------------- TERMINAL ------------------ 

Annocript will be busy checking the configuration file searching for errors and checking for correct paths and URLs. 
If all is good you will be asked something depending on what you chose in the config file and if the configuration is 
correct. 

                ---------------- TERMINAL ------------------ 
                >>> Database settings: 
                ...DONE! 
                ---------------- TERMINAL ------------------ 

Both the UniRef and the TrEMBL databases can be chosen to be downloaded. You could also install 
both (one per time) and decide later what to use during the execution. 

                ---------------- TERMINAL ------------------ 
                You chose to create the database. Do you want to create with UniRef (u) or TrEMBL (t)? (u or t) r 
                kbVersion: 2015_01 
                Folder for DB will be: /home/francesco/ann_works/data/2014_01 and database name: uniref_2015_01 
                
                The Conserved Domains database: Cdd 
                During this process Annocript will download a set of files with big dimension. These data will be used 
                to create a large  database with all the data. Do you want that Annocript erase them? (y or n). 
                Press x to exit Annocript. 
                 n 
                
                >>> Programs settings: 
                Sequences found in your session folder... 
                Annocript detected that your sequences are dna 
                >>> Setting final statistics and outputs 
                ..OK! WELL DONE! 
                Your configuration has been saved! A log file will be created in 
                /home/francesco/ann_works/jobs/new_database/log/annocript_exec_Mon_Sep_16_09:46:47_2013.log. 
                The overall computation usually depends from the number of sequences and the speed of the machine 
                you are using. 
                You will find all the results in the folder you chose here. 
                
                Please let us know if some problem occurs during the computation or you don't  retrieve what you are 
                expecting from Annocript. 
                
                ######################################### 
                # ANNOCRIPT STARTED IN BACKGROUND! BYE! # 
                ######################################### 
                ---------------- TERMINAL ------------------ 

If something goes wrong *Annocript* alerts you in the terminal output and in log file. 
Your DB will be created in few hours depending by the speed of connection and you configuration. 
You can have a look at the log file in the log folder. 

NOTE: We chose to call the database referring to the version of the UNIPROT database. 
Actually it is out the version of january 2015. We called this database uniref_2015_01. 
The date is taken automatically from a corresponding release_note file. 




##Starting Programs execution 


If you created the proteins database already and/or you want to run only analysis on your sequences then open the 
configuration file (config_user.txt) in the Annocript folder and set parameters as follows: 

                ---------------- TERMINAL ------------------ 
                francesco@Z800:~/Annocript$ nano config_user.txt 
                ---------------- TERMINAL ------------------ 

Set the name of your fasta sequences file and choose to execute only programs: 

                ---------------- TERMINAL ------------------ 
                fastaSeqs = sequences.fasta 
                
                doDbCreation = NO 
                doExecutePrograms = YES 
                doBuildOutput = NO 
                extractStatistics = NO 
                ---------------- TERMINAL ------------------ 

Choose the programs you want to execute (here we chose all): 

                ---------------- TERMINAL ------------------ 
                #What programs to execute 
                doBlastx = YES 
                doRpstblastn = YES 
                doBlastn = YES 
                doPortrait = YES 
                doDna2Pep = YES 
                ---------------- TERMINAL ------------------ 

In the configuration file for user you can set the parameters for the programs.

Put the FASTA file in the working directory (e.g. ann_works). (It will be shifted by Annocript in your jobs folder) 
Put the config_user.txt file in your working directory. 

                ---------------- TERMINAL ------------------ 
                francesco@Z800:~/Annocript$ cp config_user.txt ../ann_works/ 
                ---------------- TERMINAL ------------------ 

On command line go to your chose folder and run Annocript 

                ---------------- TERMINAL ------------------ 
                francesco@Z800:~/$ cd ../ann_works/ 
                francesco@Z800:~/ann_works/$ ../Annocript/annocript.pl config_user.txt 
                ---------------- TERMINAL ------------------ 

IMPORTANT: if you already ran some program, a second run will destroy the previous results. If some output is already 
          present, Annocript will ask if you are sure to delete. 

Annocript starts asking about the session to use. 

                ---------------- TERMINAL ------------------ 
                1. New session - 2. Use Previous - 3. Exit 
                Your choose:  1 
                Choose a name for your session: pseudonitzschia    
                /home/francesco/myAnnocript/jobs/pseudonitzschia created! 
                Copying the given config file in /home/francesco/test/jobs/pseudonitzschia 
                
                ##################################################### 
                # CHECKING THE INTEGRITY OF THE CONFIGURATION FILES # 
                ##################################################### 
                
                They are good. Annocript can proceed! 
                ################################### 
                # RUNNING ANNOCRIPT CONFIGURATION # 
                ################################### 

                Some programs will be run. In the following are listed 
                The final table will be built 
                A web page and file with statistics will be created 
                A blastx of your sequences will be done. 
                A rpstblastn of your sequences will be done. 
                A blastn of your sequences will be done. 
                Do you agree with it?(y or n) y
                
                Ok!.. let's go! 
                
                Checking : DB_Creation configuration...Checking the links... this can take a while.. please wait..... 
                ..DONE! 
                ---------------- TERMINAL ------------------ 

Since you already created the database, Annocript will ask to you what to use: 

                ---------------- TERMINAL ------------------ 
                DB Creation Not Active: Database will not be downloaded and installed. 
                You chose to use an existent database, are you sure?(y or n) y 
                
                Choose the database version: 
		2013_08  2013_09  2014_01
                ---------------- TERMINAL ------------------ 


Annocript will check if UniRef or TrEMBL databases are present and ask to you to choose one. If one only is 
present, it will be used automatically. 

                ---------------- TERMINAL ------------------ 
                Uniref database is the only present. Using it... 
                Conserved Domains database: Cdd 
                >>> Programs settings: 
                Sequences found in your session folder... 
                Annocript detected that your sequences are dna 
                >>> Setting final statistics and outputs 
                ..OK! WELL DONE! 
                Your configuration has been saved! A log file will be created in 
                /home/francesco/ann_works/jobs/pseudonitzschia/log/annocript_exec_Mon_Sep_16_09:46:47_2013.log. 
                The overall computation usually depends from the number of sequences and the speed of the machine 
                you are using. 
                You will find all the results in the folder you chose here. 
                Please let us know if some problem occurs during the computation or you don't  retrieve what you are 
                expecting from Annocript. 
                
                ######################################### 
                #STARTING ANNOCRIPT IN BACKGROUND! BYE! # 
                ######################################### 
                ---------------- TERMINAL ------------------ 


Once you did it, Annocript will start with the execution of programs that you chose in the config_user.txt file. 
If something goes wrong, Annocript will alert you in the log file. 


###Check of the results 

All BLAST output will be stored in the folder OUTPUT of your session folder. They will be checked by *Annocript* in this way: if the number of output sequences is different from the total number of sequences in your query file, than there's a problem and you will be informed. 


TO REMEMBER: You can run the programs provided by Annocript (BLAST, PORTRAIT and DNA2PEP ) also in a unique way. 
You can think to run only BLASTX and no other stuff simply because you need the output of the BLASTX. Or you can 
have the output of BLASTX for your sequences and you want only to run the other programs. 



##Print informations 


If you already get the results of the programs you can start the module to get only the statistics. 
To do this you have to inform Annocript that you want to run the modules for print out only. 
Open the configuration file in the Annocript folder and set parameters as follows: 

                ---------------- TERMINAL ------------------ 
                doDbCreation = NO 
                doExecutePrograms = NO 
                doBuildOutput = YES 
                extractStatistics = YES 
                ---------------- TERMINAL ------------------ 

(Here we do not specify to write the FASTA file name because we are considering that you already ran the 
ExecutePrograms and DBCreation modules and you are simply using the same configuration file.) 

Then run Annocript! 
                
                ---------------- TERMINAL ------------------ 
                francesco@Z800:~/ann_works ../Annocript/annocript.pl config_user.txt 
                
                1. New session - 2. Use Previous - 3. Exit 
                Your choose:  2 
                eraseme  parenysensis 
                
                Choose your folder: eraseme 
                Copying the given config file in /home/francesco/test/jobs/eraseme 
                ##################################################### 
                # CHECKING THE INTEGRITY OF THE CONFIGURATION FILES # 
                ##################################################### 
                
                They are good. Annocript can proceed! 
                ################################### 
                # RUNNING ANNOCRIPT CONFIGURATION # 
                ################################### 

                The final table will be built 
                A web page and file with statistics will be created 
                Do you agree with the following configuration?(y or n) y 
                Ok!.. let's go! 

                >>> Database settings (if you chose to create the database, it will need a while): 
                DB Creation Not Active: Database will not be downloaded and installed. 
                You chose to use some existent db files, are you sure?(y or n) y 
                2013_07 
                Choose your folder: 2013_07 
                Uniref and Uniprot databases are both there. What do you want to use? UNIREF (r) or UNIPROT (p)? 
                (r or p) r 
                The Conserved Domains database used will be: Cdd 
                
                >>> Programs settings: 
                Path to control for fasta: /home/francesco/test/trial_transcriptome.fasta 
                Sequences found in your session folder... 
                Annocript detected that your sequences are dna 
                
                >>> Setting final statistics and outputs 
                ..OK! WELL DONE! 
                
                Your configuration has been saved! A log file will be created in 
                /home/francesco/ann_works/jobs/eraseme/log/ annocript_exec_Mon_Sep_30_17:38:05_2013.log. 
                The overall computation usually depends from the number of sequences and the speed of the machine 
                you are using. 
                You will find all the results in the folder you chose here. 
                
                Please let us know if some problem occurs during the computation or you don't  retrieve what you are 
                expecting from Annocript. 
                ########################################## 
                # STARTING ANNOCRIPT IN BACKGROUND! BYE! # 
                ########################################## 
                ---------------- TERMINAL ------------------ 


The output and statistics will be all stored inside the 'stats' and 'output' folders. 


##HINTS & TIPS


Here you'll find some tricks to obtain that Annocript does something that you want but actually its normal functioning 
does not permits it.

###How to blast against specific organisms

To blast against a specific organism or a set you might use the option *blastedOrganism*. This variable can have one of two values: it can be a
file with a list of organisms or the string 'all'. The latter case represents the usual running of *Annocript* indexing all the database and blasting
against all the UniProt organisms. The first case permits to blast against a list of selected organisms inside the Trembl and SwissProt databases.

If you want that *Annocript* align your sequences against selected organisms please do the following:
- write the name of the organisms in a text file with a name which will be specific of the analysis you are doing (i.e. your_organisms). You may want to search the
available organisms by running *Annocript* and selecting '4. Do some utils' -> '4. Search organism'. With this option you can read all the available
organism's list or search a particular species.
- put the text file inside the working folder 'ann_works'
- open the configuration file and set 'blastedOrganism = your_organisms' and  doDbCreation = NO
- run Annocript

Hint: to get automatically a list of organisms from Uniprot.org you can browse the taxonomy at http://www.uniprot.org/taxonomy and build queries such as http://www.uniprot.org/taxonomy/?query=ancestor%3a33090 (ancestor:"Viridiplantae [33090]"). Then you download the tab file and pick only the column with scientific names of species related with viridiplantae.

###How to run Annocript with already downloaded db files

Since Annocript is downloading from external sources the database files you could experience some problems that are
independent from Annocript code. Sometimes a good connection is not available and permission on Uniprot servers can change
making Annocript blocking on downloads or uncompressions. 
For this reason we add simple istructions here on how to run Annocript with already downloaded db files. 
To perform this task you should simply:

- access the file config_annocript.txt and download all the files that are given as links to external source. You may want to use Linux wget or other ways.
- uncompress the files

The release date file (different for Uniref and Uniprot_kb) is needed by Annocript to create a folder to use for all the data (usually YEAR_MONTH. For example 2014_02).
This folder is created by *Annocript* in the ann_works/data folder. If it has not been created yet by *Annocript*, you must create it by reading the release date in uniref90.release_note (for Uniref90 db) or reldate.txt (for Uniprot_kb db).
Then the release date file should be inserted in the parent data folder.
The final situation should be that you have:
ann_works/data/reldate.txt present
and
ann_works/data/2014_02 containing all the uncompressed files.

Then you can start Annocript's DB_CREATION and ignore all the warnings you will receive.



###Running your own output files

You may want to run by yourself the BLAST programs and then use them with *Annocript*. It is really straightforward to do simply by putting all the files in the apposite folders of the programs.

In order to do this operation follow these steps:

1. Create your output files by using BLAST programs (without *Annocript*)

2. You need then to change the name of the output. *Annocript* is coded in such a way it can read multiple output (if you made multiple blast runs for the same transcriptome). Hence, if you made a single run you can use '1' only. Names to be given are the following:

	- blastxSPOut1, for blastx (or blastp) against SwissProt output;
	- blastxUnirefOut1, for blastx (or blastp) against UniRef output (whatever is the UniRef version);
	- blastxTROut1, for for blastx (or blastp) against TrEMBL output;
	- blastnOut1, for blastn (tblastn) output against SILVA db;
	- rpstblastnOut1, for rpsBLAST (rpstblastn) output against CDD database.

	If you have outputs from different runs and you want to join them please put a number that indicates the successive outputs (i.e. blastxSPOut1, blastxSPOut2, blastxSPOut3). 

3. Put the files just modified in the corresponding output folder of your session folder (BLASTX, RPSTBLASTN, BLASTN). If you do not know what is the session folder take a look at [how-annocript-manages-folders-and-files] (https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/OUTPUT.md#how-annocript-manages-folders-and-files). If you do not have already a session you can run once *Annocript* for some other purposes: run other programs, generation of a database or you can simply choose to run nothing just to have a session.

4. Open the config_user.txt files and modify the execution in such a way the programs that you do not want to execute must have the NO value. For example, if you have output for BLASTX against SwissProt you must have doBlastxSP=NO.
5. Run *Annocript*!



###Using different Conserved Domains Databases

If you want to use another database different from the Cdd one, you can choose one at 
(ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/). To use it you have to:
change the link cdDBLink  in the Annocript configuration file to use one of the databases present 
(Pfam, Cog, Kog, Cdd, Smart, Prk, Tigr). (The extension .rps is of the file that Annocript will search to check if the 
indexed database already exists).
Run the DB creation of Annocript.

Suppose now that you already installed DB with Annocript with Smart CD db and you want to use also Pfam. What can 
you do?

You can do it but the solution is little manual. 
You should download the DB manually with wget or some other FTP program
decompress it and put all the files inside the DATA folder where you are interested it is;
Then you can normally start Annocript: it will ask to you what cd database to use.



###Getting the output from a filtered fasta file
Suppose you ran *Annocript* for your organism but you want to investigate more a set of those transcripts. You may want to get the annotation only for those transcripts while retaining the original output from the previous run. 

You can do it!! *Annocript* will perform an identical run but on those sequences only. All the output files will be overwritten (not those from the programs). All the files in the stats folder will be overwritten. Please be careful to save all these files and remove then everything from the stats folder!!

Do it by:
- creating the filtered list of fasta sequences from the original set (filt_seq.fasta)
- putting this fasta in the session folder (the previous session you saved)
- changing the config_user.txt: change the fasta file name and changing doDbCreation = NO; doExecutePrograms = NO; doBuildOutput = YES and extractStatistics = YES.
- start again *Annocript* by loading a previous session (Option 2). Then, you type 'Y' when it asks to overwrite the configuration file in the session folder you are loading.



###Re-generation of an existing database

If you do not cancelled the compressed or uncompressed files of the database, you can later recreate it if
Annocript has been updated or the system had some problem. Simply you should put the version file in the ann_works/data
folder and run Annocript database creation again. It will find that version and will not create and download existing file.

Since Annocript is coded in such a way that it can correct an existing database, you should manually remove it.
To remove a mysql database you must first log in:

**mysql -uusername -ppassword**

then drop it:

**drop database databaseName;**

To help yourself in identifiying the name of the database you want to remove you can execute the show command:

**show databases;**s

After removal, re-start Annocript.




###Using different Uniref types of databases (100, 90, 50)

Deciding to use the UniRef database will reduce the computational time for both the construction of the database and running of blastx while the results are similar to those with uniprot_kb with only a little loose in specificity. 
The default type in *Annocript* is the UniRef90 but you may want to use another among 100 or 50. 
To use those you have simply to specify the different database links inside the annocript_config.txt file.
Variables to change are: unirefDBLink and unirefVerLink

Actually (january 2015) these links for 50 and 100 are working:
50's:
```
unirefDBLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz
unirefVerLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.release_note
100's:
unirefDBLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/uniref90.fasta.gz
unirefVerLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/uniref100.release_note
```

All the results will be given in accord with the database chosen.


###How to add parameters for blast

*Annocript* has few parameters you can give in input to blast. Actually we chosen to keep few initially and see what is the exigency of the users. If you want to add more parameters please ask in the google group. It is very simple to add if you have a little experience in coding and if you are not we can think to add this parameter.

----------------------------------------

If you get some error during the installation or the running of *Annocript* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript. 
francesco -dot- musacchia -at- szn -dot- it
