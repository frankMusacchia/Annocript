# Annocript - Installation



## Automated Annocript install 

Download the zip file on the main page of *Annocript* on GitHub or click [here](https://github.com/frankMusacchia/Annocript/archive/master.zip). Then unzip the file. One folder will come out 'Annocript-master' that we suggest to rename: 

              ---------------- TERMINAL ------------------  
              frank@compaq2:~/francesco$ unzip Annocript-master.zip 
              frank@compaq2:~/francesco$ mv Annocript-master Annocript 
              ---------------- TERMINAL ------------------ 

Create a directory 'bin' for all the softwares: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco$ mkdir bin 
              ---------------- TERMINAL ------------------ 


To install Annocript we will create a folder to save all the results produced during the analysis and a folder for the programs to install. Then we will run the complete_install.pl script from the Annocript folder. 

NB: If you don't have sudo permissions you will need to manually install MySQL, PERL, PYTHON, BioPERL and R. You must follow the instructions below (Manual installation). After you followed the instructions at **Need sudo permissions**, come back here. Otherwise, if you have sudo permissions, you can continue:

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ mkdir ann_works 
	      frank@compaq2:~$ mkdir bin 
              frank@compaq2:~$ cd Annocript 
              frank@compaq2:~$ perl complete_install.pl 
              ---------------- TERMINAL ------------------
              
This will ask to you the full path to the working folder and to the bin folder you just created: 

              ---------------- TERMINAL ------------------ 
              	This script will prepare Annocript to work in your directory. 
              		  
             	Are you sure you want to enjoy this software?(y or n) 
                y 
              	Write an existing complete path where you want to play with Annocript (/home/username/ann_works): /home/francesco/ann_works

		Write an existing complete path where you want to install the programs (/home/username/bin/): /home/francesco/bin
				   
		INSTALLING ....
		....
		....
		##############################PROGRAM LINKS #################################
		Please use the following links for the programs into the configuration file:
		blast_path = /home/francesco/bin/ncbi-blast-2.7.1+/bin/
		lncRNAPredProgPath = /home/francesco/bin/CPC2-beta/bin/CPC2.py
		dna2pepPath = /home/francesco/bin/dna2pep-1.1/dna2pep.py
		##############################PROGRAM LINKS #################################
             	 Done! Now you can run Annocript from /home/francesco/ann_works! 
              ---------------- TERMINAL ------------------ 

Notice: The complete_install.pl script gives the paths for the programs that you have to add in the configuration file later. Please copy and paste them somewhere.


Now create a user account for MySQL database.


## Instruction on how to create a user account with MySQL

If you do not have experience or a system administrator to create a mysql account but you have been able to install it, you could simply execute MySQL from the command line and create an use with a password to use later in *Annocript*.

Open a terminal and access mysql providing the password:

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ sudo mysql -uroot -p
              ---------------- TERMINAL ------------------ 

You will access the mysql terminal where you can simply add a new user:

              ---------------- TERMINAL ------------------ 
              mysql> create user 'annocript'@'localhost' identified by 'password';
              ---------------- TERMINAL ------------------ 
              
Then you have to provide the user the access to the information:

              ---------------- TERMINAL ------------------ 
              mysql> GRANT ALL PRIVILEGES ON * . * TO 'annocript'@'localhost';
              ---------------- TERMINAL ------------------ 

Finally reload the privileges to finalize it:

              ---------------- TERMINAL ------------------ 
              mysql> FLUSH PRIVILEGES;
              mysql> quit
              ---------------- TERMINAL ------------------ 

Please copy your user and password somewhere because you will need them for setting *Annocript*.


If you have any problems during the installation please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript 

To see how to use the first time Annocript go to [TUTORIAL](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/TUTORIAL.md) page




## Manual installation



If you could not install the programs with the script complete_install.pl or you want to install them manually please follow the instructions below:

**Need sudo permissions**

If you don't have sudo permissions you can still install automatically BLAST, CPC2 and dna2pep but you need to open the file Annocript/DL_PROGRAMS/install_all.sh and comment the section delimited by "########NEEDS SUDO PERMISSIONS". 

Before to use complete_install.pl Annocript needs the following packages. All of them may be safely downloaded using the Ubuntu Software Center (or other similar tools in different Linux releases). If you use this kind of tools you can simply search each package with its name and install it.
The supplementary Perl MODULES may be not included in the default installation, thus install them manually.
 
  - MySQL database server (tested with ver5.5) (http://www.mysql.com/)
  - MySQL database client with an account (tested with ver5.5) (http://www.mysql.com/)
  - Perl (tested with ver5.10) (http://www.perl.org/get.html) ;
  - BioPerl: (tested with ver1.6) (http://www.bioperl.org/wiki/Getting_BioPerl); 
  - Python: (tested with ver2.7.3) (http://www.python.org/download/); 
  - R: (tested with ver3.2.5) (http://www.r-project.org/). 
    
  - Supplementary Perl MODULES: readkey.pm, Parallel/ForkManager.pm, Switch.pm
    To install them on a Debian derivate linux distribution use from a terminal: 
    - sudo apt-get install libterm-readkey-perl 
    - sudo apt-get install libparallel-forkmanager-perl 
    
    You may want to use the following Linux commands to check if you already have these perl modules:
	- perl -MParallel::ForkManager -e "print \"Module installed.\\n\";"
	- perl -MTerm::ReadKey -e "print \"Module installed.\\n\";"

**Do not need sudo permissions**

BLAST, CPC2 and DNA2PEP can be installed automatically by the script complete_install.pl. But if for some reason you can't. Follow here. 
The tested versions are in bracket. Different versions may not work with this *Annocript* version.
Download them at: 
- blast+ : ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ (tested with ver2.7.1+)
- CPC2: http://cpc2.cbi.pku.edu.cn/download.php (tested with ver2.0beta)
  libsvm-3.18 and biopython are needed to install CPC2. Please follow the instructions contained into the README file inside the package.
- dna2pep: (tested with ver1.1) (http://www.cbs.dtu.dk/services/VirtualRibosome/download.php)

## Instructions to manually install BLAST, CPC2 and DNA2PEP 


To install please execute the following steps. Start copying both in the bin directory just created. 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco$ cp Annocript/DL_PROGRAMS/CPC2-beta.tar.gz bin 
              frank@compaq2:~/francesco$ cp Annocript/DL_PROGRAMS/dna2pep-1.1.tar.gz bin 
              frank@compaq2:~/francesco$ cd bin 
              ---------------- TERMINAL ------------------ 

TWO THINGS TO REMEMBER!! 
  1. We provide up-to-date versions of CPC2 and dna2pep into the folder DL\_PROGRAMS.  dna2pep comes with a link to Python that is incorrect. We adjusted it in the corresponding dna2pep.py script but please remember this fact if you get an error. Look inside the DL_PROGRAMS folder to find both. We suggest to install BLAST+, CPC2 and dna2pep in a common bin directory of your file system (/home/francesco/bin). 

  2. BLAST programs and DNA2PEP don't need to be installed for Annocript as it uses full paths to access the executables. If you plan to use BLAST programs later you can choose a manual installation. Instead you will need to install CPC2 and its dependency. 

### CPC2 installation 


CPC2 makes use of biopython and libsvm. Please follow the instructions into the README file to install it. The installation is straightforward.


### dna2pep Installation 

It don't need to be installed, only uncompress in the bin folder you created. Here we will also change 
permissions of the .py file or *Annocript* will not find it. 
 
(If you don't use the version included in *Annocript*, please modify the first row of the dna2pep.py with the correct path for Python (i.e. /usr/bin/python)) 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco$ cd bin 
              frank@compaq2:~/francesco/bin$ tar -zxvf dna2pep-1.1.tar.gz 
              frank@compaq2:~/francesco$ cd dna2pep-1.1
              frank@compaq2:~/francesco/dna2pepp-1.1$ chmod 755 dna2pep.py 
              ---------------- TERMINAL ------------------ 

If you have any problems during the installation please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript 


### NCBI-BLAST Installation 
              

You can both uncompress the archive compiled for your architecture and put it in the /bin folder or follow the NCBI instructions for the installation. 

Once you downloaded the file from the NCBI website you have only to unpack the tar file in the bin folder:

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco/bin$ tar -zxvf ncbi-blast-2.2.27+-x64-linux.tar.gz
              ---------------- TERMINAL ------------------ 

In Ubuntu Linux you can also install it via terminal: 

              ---------------- TERMINAL ------------------ 
              sudo apt-get install ncbi-blast+ 
              ---------------- TERMINAL ------------------ 


If you have any problems during the installation please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript 



To see how to use the first time Annocript go to [TUTORIAL](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/TUTORIAL.md) page


------------------------------------------------

If you get some error during the installation or the running of *Annocript* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
