#FAQ Frequently Asked Questions


--------------------------
*Question*: How Annocript identifies lncRNAs?

*Answer*: *Annocript* can identify long non-coding RNAs by using the complete tabular output. It goes through all the fields and if all the following requests are satisfied, then the sequence is given the value lncRNA4Annocript=1
The sequence is non coding if:
- there is no homology result with Proteins, Domains and other non-coding RNAs;
- the longest orf length is less or equal than NCORFLength (default: 100 aminoacids);
- the length of the query is equal or more than NCSeqLength (default: 200 nucleotides);
- ProbToBeNonCoding is grater than or equal to NCThresh (defaul: 0.95).

--------------------------

Here you can find issues already solved at any step with *Annocript*. Please search here before to post any question on the forum or to mail the mantainer.
   

*Question*: I got this error during the Portrait installation:
        ```
        FATAL :: Unkwnown error in CAST program. If possible, please try installing PORTRAIT in another machine. 
        CAST dump...
        ```

*Answer*: CAST is written in ANSI C, so it needs libraries and a compiler that can happen are not installed in your system.
Please install all the libraries that are requested in the INSTALL page.


--------------------------


*Question*: I got this error during the Portrait installation:
       ```FATAL :: Unkwnown error in ANGLE program. If possible, please try installing PORTRAIT in another machine. ```

*Answer*: It is a missing library that you can install via terminal.
Please install all the libraries that are requested in the INSTALL page.


If you don't used the version of Portrait included in *Annocript* you have to convert their file to Unix format. Go to
terminal and type:

```
apt-get install dos2unix
then
dos2unix portrait-1.1.pl
```

------------------------


*Question*: I got this error in the log file:

```
Thu Nov 28 12:32:37 2013: Running DNA2PEP...
sh: 1: /home/shady/bioinfo/dna2pep-1.1/dna2pep.py: not found
ERROR: an error occurred while running /home/shady/bioinfo/dna2pep-1.1/dna2pep.py: 32512 at 
PROGRAMS_EXEC/execute_programs.pm line 864.
```

*Answer*: If *Annocript* didn't find dna2pep it can means or that it is not really there or that you don't has the correct 
permissions to access it. So control if the path in the config_user is correct and, if it is, change the permissions
of dna2pep.py:

```
chmod 755 dna2pep.py.
```

It can be possible that the Python interpreter is not present. You should install it and write in the dna2pep.py file
the correct path in the first line. You should open the dna2pep.py in the dna2pep folder and in the first line change 
that path with the path where you have python. Usually it is /usr/bin/python.
If you already installed maybe that path only is not correct.

--------------------------


*Question*: Long time is taken for the downloading and/or uncompression of big compressed files. What can I do?

*Answer*:  Depending on your system memory and processors *Annocript* could have any difficulties in downloading and 
uncompressing big dimensions files. This is the reason why we suggest to run it only on powerful machines. Anyway you 
can overcome this problem by seeing what is the file that is causing the problem and downloading it and uncompressing 
it with system commands and put the resulting uncompressed files inside the ann_works/data/YEAR_MONTH folder created by 
*Annocript* for the database. Then you just have to run again it as before. *Annocript* will check if the file has been 
already downloaded or uncompressed and will not repeat the operation.

To do this:
  - download the file you think is giving problems;
  - start the DB creation execution;
  - stop *Annocript* when it is just gone in background (it will create the folder for the data);
  - pick the version file (i.e uniref90.release_note) from the folder of the database you were just installing and put
    it in the data folder;
  - put in the db version file the ZIP files you previously downloaded;
  - run again the DB creation of *Annocript*.


-------------------------

*Question*: The library libstdc++5 is missing in Angle! What can I do?

*Answer*: Please install all the libraries that are requested in the INSTALL page.


--------------------------

*Question*: The library tcc:i386 is missing in Cast! What can I do?

*Answer*: Please install all the libraries that are requested in the INSTALL page.


--------------------------

*Question*: I had the following error:
```
ERROR [256]: an error occurred while extracting headers with fgrep from /home/francesco/ann_works/data/2014_01/uniref90.fasta to /home/francesco/ann_works/data/2014_01/db_headers.txt. 
```

What does it means?

*Answer*: it could be that you ran once *Annocript* and it stopped on the download or uncompression of some file. 
Then you ran it again and it is blocked unexpectedly. You should remove completely the folder 2014_01 with the data
from the ann_works/data folder. Then run it again...

---------------------------

*Question*: I got this error message: 
```
Creating a file with all the headers in /home/francesco/ann_works/data/2014_03/db_headers.txt (it makes the db construction faster but it will need a while)...ERROR [-1]: an error occurred while extracting headers with fgrep from /home/francesco/ann_works/data/2014_03/uniprot_sprot.fasta to /home/francesco/ann_works/data/2014_03/db_headers.txt: ? at DB_CREATION/db_creator.pm line 668
```

*Answer*: With older machine but also unpredictably, it can happen that a system command does not run well with big dimensions
files. *Annocript* is programmed to try more than once the commands but if it blocks is preferrable that you manually start it again.

-----------------------------

*Question*: I had the following error:
```
DBD::mysql::db commit failed: MySQL server has gone away at DB_CREATION/db_creator.pm line xxx,
<UNIGOFILE> line 52159208.
```
...and *Annocript* stopped.

*Answer*: This error happens when the MySQL server itself stops the execution. We studied that this happened for us when the
uncompression procedure (during DB_CREATION) are not working and taking a very long time. To solve you need to restart *Annocript*
but first is very important that you see the log file and check if in the point of forced stop there are some failed process
(like uncompression). If everything is good then you should only restart again *Annocript* included the database creation.
To do this please only shift the file of uniref/uniprot_kb (depending on what you are creating) version from the current
'year_month' folder to the 'data' folder. Then restart *Annocript* using the same session.


-----------------------------

*Question*: *Annocript* suddenly stopped while it was uncompressing a file, giving the error "Out of memory!"

*Answer*: Try to download the file manually by going in the directory ann_works/data/YEAR_MONTH/ and using the 
following code:
```
tar -zxvf yourfile.tar.gz if it is a ".tar.gx" file or 
gunzip file.gz if it is a ".zip" or ".gz" file
```

Then, start *Annocript* again. It will find the uncompressed file and will not do the operation automatically.


-----------------------------

*Question*: I cannot find anymore my fasta in the *Annocript* working directory (e.g. ann_works)! What is happened to my transcriptome?

*Answer*: Be aware that if you erase the working folder, also you transcriptome will be removed. We strongly suggest to always keep a copy of the transcriptome in some other folder.


-----------------------------

*Question*: I got this error: ```Fetch failed! HTTP response: 500 Internal Server Error [500 LWP::Protocol::MyFTP: connect: Connection timed out] at USEFUL/utilities.pm line 886```
during the download of files. Is there a problem?

*Answer*: No problems if the file has been  finally downloaded. *Annocript* does more than one trial to download the  files and if it does not succeed, will give an error and blocks completely the execution.


-----------------------------

*Question*:  I got this error:
```
Please check if [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz] is a correct URL. *Annocript* will continue...
Since there is some problem with connection, Annocript is trying different kind of download. Be sure that the links above are working and then press y to continue.
 Continue? (y or n)
```

*Answer*: This error comes when a problem with sending the request from PERL function is occurred. *Annocript* uses different PERL functions to download the files so you should check if this links are good an proceed with the running. If the connection is not working well *Annocript* will stop during the download of files.

-----------------------------

*Question*: When I was choosing the session to use, *Annocript* wrote this:
```Session astropecten2 is not valid. Do you want to remove it?(y or n) n.``` What happened?

*Answer*: *Annocript* controls everytime if the session you want to use is valid. The check consists in looking if all the folders constituting a valid session are present. They are: log, stats, data, output, gff. If something is missing *Annocript* will stop.

------------------------------


*Question*: *Annocript* stopped suddenly giving the following error (or with other file to download):

``` 
Downloading File: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
Some problems occurred while downloading ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz. Trying again...
Fetch failed! HTTP response: 404 Not Found [404 File 'idmapping_selected.tab.gz' not found] at USEFUL/utilities.pm line 1145
Could not fetch '/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz' from 'ftp.uniprot.org' at USEFUL/utilities.pm line 1145
``` 

What happened?

*Answer*: this kind of error means that something during the download was not successful. Since *Annocript* has two or more ways to access the files
on internet we may consider that something is happening on the server were we are trying to download. To be sure that the reason is this, you may 
access the link that is giving the problem with a normal browser or with **wget** from the terminal. If you can access the file in this way, then please refer to our mailing list. 
Otherwise it is a problem that is actually affecting the server. You should wait until the server is up and then you can restart *Annocript*.

------------------------------

*Question*: *Annocript* gives always this "Fetch failed" error, even do then it continues and downloads quietly... Why this happens?

*Answer*: As we said before, *Annocript* tries different ways to download a file and more than one time. It can happen that the server is busy or that the request is not
well sent and the server can not to respond in the expected way. Thus *Annocript* tries different times and functions to access and download the files.
Depending also on your internet connection this error can come or not and you can be confused when you read it on the log file. Please check always if the file is currently
being downloaded in the data/YEAR_MONTH folder before to start again *Annocript*.

------------------------------

*Question*: *Annocript* stopped with this error:

```
Timeout at /usr/share/perl5/LWP/Protocol/ftp.pm line 278 at /usr/share/perl5/LWP/UserAgent.pm line 844.
```

what does it means?

*Answer*: While it was downloading a file the server interrupted the comunication. It is a rare situation but since it can happen you should remove the partially
downloaded files and restart *Annocript*. The pipeline will download only the files that are missing.


------------------------------

*Question*: *Annocript* stopped with this error:

```
perl: warning: Falling back to the standard locale ("C").
Segmentation fault (core dumped)
Illegal division by zero at /home/francesco/bin/portrait-1.1/portrait-1.1.pl line 608.
ERROR [35584]: an error occurred while running /home/francesco/bin/portrait-1.1/portrait-1.1.pl: 35584 at PROGRAMS_EXEC/execute_programs.pm line xxx
```

what does it means?

*Answer*: This error is common with Portrait. It is a stochastic interruption of the analysis we have not been able to solve in a direct clean way.
Anyway we found a trick to let Portrait work also in this case! We saw that when group of sequences are used with this program 
instead of all together, it works. Thus we created a variable in the configuration file which you are actually using, that you can
use to tell Annocript to separate your N sequences in groups of M sequences

-------------------------------

*Question*: *Annocript* stopped with this error:

```
FATAL :: Output file /data02/francesco/ann_works/jobs/human_mitoblaes/portrait_human_mitoplast_transcriptome_5.fasta_results already exists!
```
what does it means?

*Answer*: This happens when you have a false Portrait run. Portrait creates its own temporary files in you working folder. Then if an erroneous
execution happen for any reason, the temporary files are not removed. The solution is to remove all the files containing
the word "portrait" in your working folder. Then run again Annocript. Use the configuration file to instruct Annocript
to run only the programs that failed! 

-------------------------------

*Question*: *Annocript* stopped with this error:

```
The type of UniRef you want to download is not the same as the one you downloaded in a previous execution of Annocript. Annocript can manage only one type of Uniref (50,90 or 100) for each versions of the Uniprot db.Check the unirefDBLink and unirefVerLink in CONFIGURATION/config_annocript.txt file.
Died at ../Annocript/annocript.pl line xxx.
```

*Answer*: This stop of *Annocript* arises when the user downloads a Uniref database (say Uniref50) and later he wants to install
also Uniref90. *Annocript* is designed in such a way it detect what is the type of Uniref installed and can install only one type
for each version. Thus you must modify the parameter unirefDBLink and unirefVerLink in CONFIGURATION/config_annocript.txt file 
or if you want to use Uniref90 and loose Uniref50, which you have installed, you need to remove the database uniref_YEAR_MONTH and
remove the Uniref90.fasta file contained in the YEAR_MONTH DB folder (the zipped file too, if you kept it). *Annocript* will download
only those files that are not contained in the folder. 

-------------------------------

-------------------------------

*Question*: *Annocript* stopped with this error:

```
trying URL 'http://cran.us.r-project.org/src/contrib/RColorBrewer_1.0-5.tar.gz'
Error in download.file(url, destfile, method, mode = "wb", ...) : 
  cannot open URL 'http://cran.us.r-project.org/src/contrib/RColorBrewer_1.0-5.tar.gz'
In addition: Warning message:
In download.file(url, destfile, method, mode = "wb", ...) :
  unable to connect to 'cran.us.r-project.org' on port 80.
Warning in download.packages(pkgs, destdir = tmpd, available = available,  :
  download of package ‘RColorBrewer’ failed
> #else{print("RColorBrewer already installed")}
> 
> 
> library(RColorBrewer)
Error in library(RColorBrewer) : 
  there is no package called ‘RColorBrewer’
```

*Answer*: *Annocript* was not able to install a necessary library to add colors to your plots. Please try to install it manually going in the R console and typing: install.packages('RColorBrewer'). Then restart *Annocript* with only extractStatistics = YES


---------------------------------------------------------------------------------

If you get some error during the installation or the running of *Annocript* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
