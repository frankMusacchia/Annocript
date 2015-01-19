#Annocript 0.2 - Configure and run!

##Configuration 

The first time you run *Annocript* you have to fix some important parameters in the config_user.txt file. 
These parameters will not change anymore after you set them the first time. Then you can keep that configuration file in your working folder and modify it everytime you need a different analysis. The configuration file is in the CONFIGURATION folder of *Annocript*. Copy that and open it with a text editor. I use nano:

                ---------------- TERMINAL ------------------
                francesco@compaq2:~/ann_works$ cp ../Annocript/CONFIGURATION/config_user.txt .
                francesco@compaq2:~/ann_works$ nano config_user.txt 
                ---------------- TERMINAL ------------------ 


**Fasta transcriptome** 

The first parameter is the name of your FASTA file: 

        fastaSeqs = your_transcriptome.fasta
        
IMPORTANT: The name of the file with sequences must respect these constrains:        
- the extension must be .fasta or .fa;
- the overall lenght (with  extension) must be no longer than 50 chars;
- the name has not to contain dots. Only the one which separates the extension is accepted;
- chars accepted in the name are: A-z,a-z,0-9,_,-.
- the fasta file must be into the working directory (ann_works)

If these constrains are not respected *Annocript* will die and you will have to change the sequences file names. 

**Organisms to blast with**

This parameter is useful to define which organisms to blast with. Here you can both give a file name with a list of organisms names
or select 'all', to intend that you want to use all the UniProt database. 

If you want to blast against selected organisms you must use TrEMBL database. Look [here](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/ADVANCED_USAGE.md#how-to-blast-against-specific-organisms) for more info...

blastedOrganism = all

**GO terms associations**

You can ask *Annocript* to display GO terms associated with proteins, domains or both. To leverage this behaviour use the parameter:
**goTermsAss**. You can give it the following values: 'proteins', 'domains' or 'both'.


Reach the final part of the file (FIXED PARAMETERS) and do the following: 

**Database account info** 
You have to replace here your account informations to access the MySQL database. 
        mySqlUser  //Username for MySQL in your system 
        mySqlPass  //Password for MySQL in your system 


IMPORTANT: To use Annocript you need access to a MySQL server with an account with grants to create databases
 and tables and insert/modify/delete rows. If you don't have, please contact your system administrator. 
 

**Uniprot user access** 
Here you choose the account name and password to access and download from UniProt: 

        uniprotWebUser = anonymous (you can leave this) 
        uniprotWebPass = youremail@the_domain.com 


**Programs paths** 
Here you need to add paths to the programs that depend on the folder where you chose to install. 

        #Programs Paths 
        blastPath = /home/francesco/bin/blast/bin/ 
        portraitPath = /home/francesco/bin/portrait/portrait-1.1.pl 
        dna2pepPath = /home/francesco/bin/dna2pep/dna2pep.py 

**Using GFF DB**

This variable permits the generation of a GFF database containing all the output from the programs. Use YES only if you think you will need the GFF database. GFF files will be always print.
        
        useGFFDB = NO 
        
        
Now you are ready to use *Annocript*! 

In the following is shown how to make different things. You should first create a database for proteins and then 
run the programs. If you need only to do one of the three steps you can go to ADVANCED_USAGE page. 

READ CAREFULLY: 
- Depending upon your network you may experience some problems during the automatical download and uncompression of Uniprot files or other flat files used as sequences/annotation sources. 
Some Warnings also may come out because of problems with connection on your client and db servers. If you meet some of these problems during the Annocript run we suggest follow carefully instructions on "How to run Annocript with already downloaded db files" in the HINTS AND TIPS section in ADVANCED_USAGE page.
- Database files links are contained in the annocript_config.txt file. It is inside the Annocript/CONFIGURATION folder. The default UniRef database is UniRef50 but you may want to choose 90 or 100, too. If you want, please change both the links unirefDBLink and unirefVerLink accordingly (you need only to substitute the number!).
- The default database for domains is the Cdd. It contains domains from KOG, COG, Pfam, TIGR, Prk. Some problem could be verified when you are using a machine with 32bits. Please consider to use only one single database. We suggest Pfam but you can use which you want. If you do that, please change the variable cdName4Expression in the same annocript_config.txt file.

        
##Annocript fast complete execution 


After that the fixed parameters have been set, you can start *Annocript* and its analysis. 
Put the file with your sequences inside your working directory folder (ann_works).
A trial file (trial\_transcriptome.txt) is present in the folder USEFUL of the program.
Suppose the sequence are in the file trial\_transcriptome.txt: it must be placed in your session folder (e.g. ann\_works). 

Then open the user config file: 

                ---------------- TERMINAL ------------------ 
                francesco@compaq2:~/ann_works$ nano config_user.txt 
                ---------------- TERMINAL ------------------ 

and set the following: 

                ---------------- TERMINAL ------------------ 
                doDbCreation = YES 
                doExecutePrograms = YES 
                doBuildOutput = YES 
                extractStatistics = YES 
                ---------------- TERMINAL ------------------ 

You should also set the number of threads you want to use (depending from the number of processors you have on your system) and the name of the file with transcripts. 

                ---------------- TERMINAL ------------------ 
                fastaSeqs = trial_transcriptome.fasta 
                .... 
                .... 
                #Number of threads for parallel executions (Actually used for RPSBLAST) 
                threads4Parallel = 10
                ---------------- TERMINAL ------------------ 

Then close the editor and run Annocript: 

                ---------------- TERMINAL ------------------ 
                francesco@compaq2:~/ann_works$ perl ../Annocript/annocript.pl config_user.txt 
                ---------------- TERMINAL ------------------ 

*Annocript* will ask to you to choose something. Please keep at hand this guide while you are running our software 
to know how to answer. 

                ---------------- TERMINAL ------------------ 
                1. New session - 2. Use Previous - 3. Exit 
                Your choose: 1 
                Choose Name: myProject 
                
                A folder will be created with this new name. 
                
                Copying the given config file in /home/francesco/ann_works/jobs/myProject 
                ##################################################### 
                # CHECKING THE INTEGRITY OF THE CONFIGURATION FILES # 
                ##################################################### 
                
                They are good. Annocript can proceed! 
                ################################### 
                # RUNNING ANNOCRIPT CONFIGURATION # 
                ################################### 
                
                DB file will be downloaded from UNIPROT website and DB built. 
                Some programs will be run. In the following are listed 
                A BLASTx  of your sequences will be done. 
                A rpsBLAST of your sequences will be done. 
                A BLASTn  of your sequences will be done. 
                PORTRAIT will be executed with your sequences. 
                DNA2PEP will be executed with your sequences. 
                The final table will be built 
                A web page and file with statistics will be created 
                Do you agree with it?(y or n) 
                y 
                ---------------- TERMINAL ------------------ 


(NOTICE: if you get some error during one Session of work you can start again and load the same session) 

Once you accepted the configuration *Annocript* will create folders and will control your parameters: 

                ---------------- TERMINAL ------------------ 
                /home/francesco/ann_work/jobs/myProject/output/blastx doesn't exists. Creating it... 
                /home/francesco/ann_work/jobs/myProject/output/rpstblastn doesn't exists. Creating it... 
                /home/francesco/ann_work/jobs/myProject/output/blastn doesn't exists. Creating it... 
                /home/francesco/ann_work/jobs/myProject/gff doesn't exists. Creating it... 
                
                Database settings (if you chose to create the database, it will need a while): 
                Checking the links...DONE! 
                You chose to create the database. Do you want to create with UNIREF (r) or UNIPROT (p)? (r or p) r 
                ---------------- TERMINAL ------------------ 

Once you choose the database to use *Annocript* will download the version and set the folder for the data: 

                ---------------- TERMINAL ------------------ 
                 Downloading File: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note 
                ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note downloaded successfully! 
                kbVersion: 2013_10 
                Folder for DB will be: /home/francesco/ann_works/data/2013_10 and database name: uniref_2013_10 
                /home/francesco/ann_work/data/2013_10 doesn't exists. Creating folder /home/francesco/ann_work/data/2013_10... 
                Conserved domains database: Cdd_LE 
                
                During this process Annocript will download a set of files with big dimensions. These data will be used to create a large database with all the data. Do you want that Annocript erase them? (y or n). Press x to exit Annocript. 
                n 
                ---------------- TERMINAL ------------------ 

If you don't delete the compressed file, if there will be some error, you'll not need to download them again. We suggest to do this the first time you create a database. 

Other controls will be performed on parameters of programs execution and statistics: 

                ---------------- TERMINAL ------------------ 
                Programs settings: 
                Path to control for fasta: /home/francesco/ann_works/trial_transcriptome.fasta 
                Sequences found in your working folder. Moving in session folder... 
                Annocript detected that your sequences are transcripts 
                Number of sequences: 333 
                 
                Setting final statistics and outputs 
                ..OK! WELL DONE! 
                Your configuration has been saved! A log file will be created in 
                /home/francesco/ann_works/jobs/myProject/log/annocript_exec_Thu_Nov__7_12:56:42_2013.log. 
                The overall computation usually depends from the number of sequences and the speed of the machine you 
                are using. 
                You will find all the results in the folder you chose here. 
                
                Please let us know if some problem occurs during the computation or you don't  retrieve what you are 
                expecting from Annocript. 
                ########################################## 
                # STARTING ANNOCRIPT IN BACKGROUND! BYE! # 
                ########################################## 
                ---------------- TERMINAL ------------------ 


At this point *Annocript* is running! 


NOTE: We chose to call the database referring to the version of the UniProt database. 
Actually it is out the version of june 2014. We called this database uniref_2014_06. 
The date is taken automatically from a corresponding release_note file at:

    uniprotVerLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt

for UniProt_kb and at

    unirefVerLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note

for UniRef db.


You can look at the log file to see the current state with the linux program 'more'. 

                ---------------- TERMINAL ------------------ 
                /home/francesco/ann_works/jobs/myProject/log/annocript_exec_Thu_Nov__7_12:56:42_2013.log. 
                ---------------- TERMINAL ------------------ 

You may want to know if *Annocript* is actually running. Watch this by using the 'htop' linux command. You should see the command "perl PROGRAMS\_EXEC/annocript_executor.pl path/to/yourSequences.fasta" running. To search it press the F3 key and search 'annocript'. It should appear!

You now may want to understand what is coming out from *Annocript*. Read it in the [OUTPUT](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/OUTPUT.md) page.
Or you may want to do something more specific or change parameters, then look at [ADVANCED_USAGE](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/ADVANCED_USAGE.md) page.


------------------------------------------


If you get some error during the installation or the running of Annocript please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
