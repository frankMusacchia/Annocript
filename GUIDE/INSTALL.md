# Annocript 0.2 - Installation


## Install Dependencies

To work, Annocript need the following packages. All of them may be safely downloaded using the Ubuntu Software Center (or other similar tools in different Linux releases). If you use this kind of tool you can simply search each package with its name and install it.
The supplementary Perl MODULES bay not be included in the default installation, thus install them manually.
 
  - MySQL database server (tested with ver5.5) (http://www.mysql.com/)
  - MySQL database client with an account (tested with ver5.5) (http://www.mysql.com/)
  - Perl (tested with ver5.10) (http://www.perl.org/get.html) ;
  - BioPerl: (tested with ver1.6) (http://www.bioperl.org/wiki/Getting_BioPerl); 
  - Python: (tested with ver2.7.3) (http://www.python.org/download/); 
  - R: (tested with ver3.1.1) (http://www.r-project.org/). 
    
  - Supplementary Perl MODULES: readkey.pm, Parallel/ForkManager.pm, Switch.pm
    To install them on a Debian derivate linux distribution use from a terminal: 
    - sudo apt-get install libterm-readkey-perl 
    - sudo apt-get install libparallel-forkmanager-perl 
    
    You may want to use the following Linux commands to check if you already have these perl modules:
	- perl -MParallel::ForkManager -e "print \"Module installed.\\n\";"
	- perl -MTerm::ReadKey -e "print \"Module installed.\\n\";"
  


Then you have to download and install the programs that Annocript uses: BLAST, PORTRAIT and DNA2PEP. The tested versions are in bracket. Different versions may not work with this *Annocript* version.  
Download them at: 
- blast+ : ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ (tested with ver2.2.30)
- Portrait: http://bioinformatics.cenargen.embrapa.br/portrait/ (tested with ver1.1)
  The following libraries are needed for Portrait: g++, tcc:i386, libstdc++5. Please be sure they are present or download them with the following commands:
    - sudo apt-get install g++
    - sudo apt-get install tcc:i386 (or install both tcc and libc6-dev-i386)
    - sudo apt-get install libstdc++5
   
    You may want to use the following Linux commands to check if you already have these Linux packages:
    dpkg -l |grep g++
    dpkg -l |grep tcc:i386
    dpkg -l |grep tcc:i386

- dna2pep: (tested with ver1.1) (http://www.cbs.dtu.dk/services/VirtualRibosome/download.php)

Read carefully the instructions below to install them. 

TWO THINGS TO REMEMBER!! 
  1. Since we experienced problems with Portrait on Linux based systems when analyzing very large and complex set of sequences, we inserted a working version inside the folder DL\_PROGRAMS. This does not change the algorithm, it simply affects the way the sequences are handled. Therefore you can avoid to download PORTRAIT from the program link. Moreover, dna2pep comes with a link to Python that is incorrect. We adjusted it in the corresponding dna2pep.py script but please remember this fact if you get an error. Look inside the DL_PROGRAMS folder to find both. We suggest to install BLAST+, Portrait and dna2pep in a common bin directory of your file system (/home/francesco/bin). 

  2. BLAST programs and DNA2PEP don't need to be installed for Annocript as it uses full paths to access the executables. If you plan to use BLAST programs later you can choose a manual installation. Instead you will need to install PORTRAIT and its dependencies. 

## Instruction on how to create a user account with MySQL

If you do not have experience or a system administrator to create a mysql account but you have been able to install it, you could simply execute MySQL from the command line and create an use with a password to use later in *Annocript*.

Open a terminal and write the following:

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ sudo mysql 
              ---------------- TERMINAL ------------------ 

You will access the mysql terminal where you can simply add a new user:

              ---------------- TERMINAL ------------------ 
              mysql> create user 'frank'@'localhost' identified by 'password';
              ---------------- TERMINAL ------------------ 
              
Then you have to provide the user the access to the information:

              ---------------- TERMINAL ------------------ 
              mysql> GRANT ALL PRIVILEGES ON * . * TO 'frank'@'localhost';
              ---------------- TERMINAL ------------------ 

Finally reload the privileges to finalize it:

              ---------------- TERMINAL ------------------ 
              mysql> FLUSH PRIVILEGES;
              mysql> quit
              ---------------- TERMINAL ------------------ 

Please write down this setting to use later in Annocript.


## Instructions to install BLAST, PORTRAIT and DNA2PEP

Download the file on the main page of *Annocript* on GitHub or click [here](https://github.com/frankMusacchia/Annocript/archive/master.zip). Then go to the folder where you downloaded it, suppose 
"/Download", and unzip the file. One folder will come out 'Annocript-master' that we suggest to rename: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/Downloads$ cp Annocript-master.zip ../Annocript-master.zip 
              frank@compaq2:~/Downloads$ cd .. 
              frank@compaq2:~/francesco$ unzip Annocript.zip 
              frank@compaq2:~/Downloads$ mv Annocript-master Annocript 
              ---------------- TERMINAL ------------------ 

Create a directory 'bin' for all the softwares: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco$ mkdir bin 
              ---------------- TERMINAL ------------------ 

Now you will enter the *Annocript* directory to configure and install Portrait and DNA2PEP: 

To install please execute the following steps. Start copying both in the bin directory just created. 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco$ cp Annocript/DL_PROGRAMS/portrait-1.1.tar.gz bin 
              frank@compaq2:~/francesco$ cp Annocript/DL_PROGRAMS/dna2pep-1.1.tar.gz bin 
              frank@compaq2:~/francesco$ cd bin 
              ---------------- TERMINAL ------------------ 


### Portrait installation


Portrait makes use of 3 different programs to run: CAST, ANGLE, LIB_SVM. Each of the packages needs to be present and 
working. Cast and Angle are simply two folders with executable files. You need only to uncompress and place them in the 
Portrait folder. 

Unpack portrait.tar.gz:  

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco/bin$ tar -zxvf portrait-1.1.tar.gz 
              ---------------- TERMINAL ------------------ 

You will find the three programs in the Portrait folder. So, go there: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco/bin$ cd portrait-1.1 
              frank@compaq2:~/francesco/bin/portrait-1.1$ 
              ---------------- TERMINAL ------------------ 

Of the three packages in Portrait's folder only libsvm must be compiled. Go in the corresponding folder and compile it: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco/bin/portrait-1.1$ cd libsvm-2.84 
              frank@compaq2:~/francesco/bin/portrait-1.1/libsvm-2.84$ make 
              ---------------- TERMINAL ------------------ 

Then go in the portrait folder and type: 

              ---------------- TERMINAL ------------------ 
              frank@compaq2:~/francesco/bin/portrait-1.1/libsvm-2.84$ cd .. 
              frank@compaq2:~/francesco/bin/portrait-1.1$ perl portrait-1.1.pl 
              ---------------- TERMINAL ------------------ 
              
The script for the configuration will run. This is how it will appear: 

              ---------------- TERMINAL ------------------ 
              Welcome! Now starting configuration process. 
              Directory where PORTRAIT is located (e.g. /home/user/program/): /home/francesco/bin/portrait-1.1/      
              Directory where LIBSVM is installed (e.g. /home/user/libsvm-2.84/): /home/francesco/bin/portrait-1.1/libsvm-2.84 
              Directory where CAST is installed (e.g. /home/user/cast-linux/): /home/francesco/bin/portrait-1.1/cast-linux  
              Directory where ANGLE is installed (e.g. /home/user/angle/): /home/francesco/bin/portrait-1.1/angle 
              Running ANGLE 32-bit version. 
              Config process successful! Run again to start working. 
              ---------------- TERMINAL ------------------ 
              
If you have any problems during the installation please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript

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


## Install Annocript


To install Annocript we will create a folder where to save all the files and folder produced during the analysis. Then we will run the install.pl script from 
the Annocript folder. 


              ---------------- TERMINAL ------------------ 
              frank@compaq2:~$ mkdir ann_works 
              frank@compaq2:~$ cd Annocript 
              frank@compaq2:~$ perl install.pl 
              ---------------- TERMINAL ------------------
              
This will ask to you the full path to the working folder you just created: 

              ---------------- TERMINAL ------------------ 
              This script will prepare Annocript to work in your directory. 
                
              Are you sure you want to enjoy this software?(y or n) 
                y 
              Write an existing complete path where you want to play with Annocript: /home/francesco/ann_works/ 
                
              Done! Now you can start Annocript from this folder! 
              ---------------- TERMINAL ------------------ 

Notice: The install.pl script simply creates a file (folders.txt) where are written the folders that Annocript will use. 
If the folders.txt file is not there Annocript cannot work. 

To see how to use the first time Annocript go to [TUTORIAL](https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/TUTORIAL.md) page


------------------------------------------------

If you get some error during the installation or the running of *Annocript* please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
