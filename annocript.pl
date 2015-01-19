#!/usr/bin/perl

#Annocript - A complete tool for transcriptomes annotation
#Copyright (C) <2014>  <Francesco Musacchia>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    
=head1 NAME
Annocript - A complete tool for transcriptomes annotation
   

=head1 SYNOPSIS

Annocript 0.2 is a tool for the annotation of transcriptomes.

Please read carefully the user Readme and User Guide before to start executing it.
 
You must give in input the configuration file to run Annocript:

annocript.pl config_user.txt

=cut


=head1 OPTIONS

=over 8

=item B<--help>
Shows the brief help informations

=item B<--version>
Shows the version

=back

=cut


=head1 EXAMPLES

  Write the config.txt file and give it in input to this script. In the config file you have to specify the locations of the data
  that you want to manage.

  The following is an example of call to this script:

 annocript.pl config.txt

=cut


=head1 DESCRIPTION
  
  Here you have the descriptions of the elements that you have to write on the file config.txt
 
             
   
=cut


=head1 AUTHORS

 Francesco Musacchia
 Swaraj Basu
 Remo Sanges

 
 --
 http://www.szn.it

 $Id: Annocript,v 0.2 2012/14/11 11:00:12 Exp $
=cut

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Cwd;#To change work directory
 
#USE THE GENERAL UTILITY MODULES
use File::Copy;#To manage file
use FindBin;
use lib $FindBin::Bin;
use Term::ReadKey;
use Getopt::Long;#To control the input
use Pod::Usage;#Used to write an usage
use LWP::Simple;#Functions for FTP management
use File::Path qw(make_path remove_tree);#To remove directories
use IO::Handle;#To immediately print with autoflush 
 
#Annocript libraries
use DB_CREATION::db_creator qw( execute_db_creator);
use PROGRAMS_EXEC::execute_programs qw( execute_programs  );
use GFF3_AND_OUTPUT::gff3_manager qw(execute_gff_manager);
use USEFUL::utilities qw(checkLink correct_type db_present extract_name 
                  detect_fasta_type clean_cd_name print_array 
                  num_processors append_file_2_file check_presence my_head ram_memory
                  is_folder_empty save_hash testFTP check_FTP_diff_sources indexed_db_present getDBVersion);
use USEFUL::utils qw(execute_utils);

#EXTRACT PARAMETERS FROM THE CONFIG FILE
my $configHash;#This hash substain all the life of Annocript
my $version = '0.2.29';
my $dbTablesVersion = '0.2';
my $configUser = "config_user.txt";
my $configAnnocript = "CONFIGURATION/config_annocript.txt";
my $variablesFile = "CONFIGURATION/variables.txt";
my $dataFolder = "data";
my $logFolder = "log";
my $sessionFolder = "jobs";
my $outFolder ="output";
my $gffFolder = "gff";
my $statsFolder = "stats";
my $tempTableFolder = 'tables';#Folder to keep the tables while they are being generated

my $usefulFolder = 'USEFUL';#A folder with useful scripts and data
my $annocriptFolder;#Folder with scripts
my $workingFolder;#Folder where Annocript has to work
my $foldersFile = "folders.txt";
my $logFile = "";
my $timesFile = "";

#Variables present also in utils.pl (IF YOU CHANGE IT HERE CHANGE ALSO IN UTILS.PL!!)
my $jobsFolder = "jobs";#Folder Users 
my $mainDataFolder = "data";#Folder inside Annocript with all data needed more than one time

my $utilsScript = "USEFUL/utils.pl";

#We need these variables to check if the programs are present
my @blastProgramsInUse = qw( blastx
                             blastp
                             blastn
                             tblastn
                             rpstblastn
                             rpsblast
                             makeblastdb
                             blastdbcmd
                            );
                            
#We need these variables to check if the files of the indexed cd database are present                            
my @cdFilesExts = qw(.rps .psq .psi .psd .pin .phr .loo .aux);      

#Admitted database in use for conserved domains    
my @cdDBAdmitted  = qw(cd pfam kog cog tigr smart prk);

######################################################### MAIN ########################################################
#  Parse command line arguments.  
parse_command_line_args();


#LOG FILE CREATION
#Here we take the current time to get a log name different at each computation  
my $time = scalar(localtime);
$time =~ s/ /_/g;
$time =~ s/:/-/g;
$logFile = "annocript_exec_$time.log";
$timesFile = "annocript_times_$time.log";

#open the log file
my $logFHandle;
open($logFHandle, ">$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";

STDOUT->autoflush(1);#This makes STDOUT hot in the sense that everything will be print immediately

annoPrint ("***********************************************************************\n");
annoPrint ("This program comes with ABSOLUTELY NO WARRANTY; for details http://opensource.org/licenses/GPL-3.0.
This is free software, and you are welcome to redistribute it
under certain conditions; http://opensource.org/licenses/GPL-3.0 for details.\n");
annoPrint ("***********************************************************************\n");
annoPrint ("  /&&&&&&                                              /&&           /&&           /&&&&&&     /&&&&&&   \n");
annoPrint (" /&&__  &&                                            |__/          | &&          /&&&_  &&   /&&__  &&  \n");
annoPrint ("| &&  | &&/&&&&&&& /&&&&&&&   /&&&&&&  /&&&&&&& /&&&&&& /&& /&&&&&& /&&&&&&      | &&&&| &&  |__/  | &&  \n");
annoPrint ("| &&&&&&&| &&__  &| &&__  && /&&__  &&/&&_____//&&__  &| &&/&&__  &|_  &&_/      | && && &&    /&&&&&&/  \n");
annoPrint ("| &&__  &| &&  | &| &&  | &|  &&  | &| &&     | &&  |__| &| &&  | && | &&        | &&| &&&&   /&&____/   \n");
annoPrint ("| &&  | &| &&  | &| &&  | &|  &&  | &| &&     | &&     | &| &&  | && | && /&&    | && | &&&  | &&        \n");
annoPrint ("| &&  | &| &&  | &| &&  | &|   &&&&&&|  &&&&&&| &&     | &| &&&&&&&/ |  &&&&/    |  &&&&&&//&| &&&&&&&&\n");
annoPrint ("|__/  |__|__/  |__|__/  |__/ |______/ |_______|__/     |__| &&____/   |___/      \______/__/ |________|\n"); 
annoPrint ("                                                          | &&\n");                        
annoPrint ("                                                          | &&\n");                                           
annoPrint ("                                                          |__/\n");  
annoPrint ("\nAnnocript is a tool created by the Bioinformatic Lab of SZN of Naples.");# With this software you can execute a series of commands");
annoPrint (" \nGiven a fasta file with a transcriptome, it will annotate your sequences and");
annoPrint (" it will separate putative coding and long non-coding RNA sequences.");
annoPrint ("\n***********************************************************************\n\n\n");


#Creates some folders that Annocript will use for the session
createBasicFolders();

#Choose what session or if one new
my $newSession = showOptions();

#We check if the config file is already present in the given folder. If this happens, the user can decide the overwriting
find_config_file();

#The variable will take the complete path as name
$configUser = $sessionFolder."/".$ARGV[0];

#Here we change the working directory
chdir $annocriptFolder;

#prints a sentence in a square of asteriscs 
nice_printing("CHECKING THE INTEGRITY OF THE CONFIGURATION FILES");  

#Checks if the configuration files are present and if they have all the variables correctly written
checkConfigVariables($configUser,$configAnnocript,$variablesFile);

#Open the config file and builds the hashes for programs
configFile2Hash($configUser); 
annoPrint ("They are good. Annocript can proceed!\n");
  
#Open the config file and builds the hashes for programs
configFile2Hash($configAnnocript);

#This routine simply add some other elements to the configuration Hash
nice_printing("RUNNING ANNOCRIPT CONFIGURATION");

#Checks folder and files if they are well written
check_Annocript_integrity();

#Complete the configuration by controlling every variable inserted in the config file and  
#creating folders inside the session one
completeConfiguration();
  
 
#Saves the hash in a file to be used later
my $hashName = "configHash";
my $configHashPath = $dataFolder."/".$hashName;
save_hash($configHashPath,$configHash, $hashName);  

annoPrint ("Your configuration has been saved! A log file will be created in $logFolder/$logFile. \nThe overall computation usually"
      ." depends from the number of sequences and the speed of the machine you are using. \nYou may want to check Annocript step-by-step "
      ."with the following command: more /home/francesco/ann_works/jobs/trial/log/annocript_exec_Tue_Dec_23_09-36-50_2014.log.\n"
      ."\nYou will find all the results in the folder you chose here. "
      ."\n\nPlease let us know if any problem occurs during the computation or you don't get from Annocript the expected result.\n"
      ."You may want to use the forum at: https://groups.google.com/forum/#!forum/annocript.");

nice_printing("STARTING ANNOCRIPT IN BACKGROUND! BYE!");

close ( $logFHandle );


#print "moving the file now...";
my $newLogFile =  $logFolder."/".$logFile;#Log of Annocript goes in the log folder
move($workingFolder.'/'.$logFile, $newLogFile) or annoDie(" Annocript will not move $workingFolder/$logFile to $newLogFile\n");
$logFile = $newLogFile;
$timesFile = $logFolder."/".$timesFile;#File with times goes inside the log folder


#We also remove all the other remaining log file in the folder created by user's interruptions
system("rm -f $workingFolder/annocript_exec_*"); 
    
            
#This perl script will be run in background
my $runAnnocriptCmd = "nohup perl ".$configHash->{'ProgExecFolder'} ."/annocript_executor.pl $configUser $configAnnocript $timesFile $configHashPath >> $logFile &";
( system($runAnnocriptCmd) ) == 0 
  or die("Unable to start annocript_executor!\n");

###############################################################CHECKS AND CONTROLS###########################################


=head2 find_config_file

 Title   : find_config_file
 Usage   : find_config_file()
 Function: check if a config file is present in the folder. If yes, it can overwrite it or not
            
 
 Returns : nothing

=cut
sub find_config_file{
  annoPrint ( "Copying the given config file in $sessionFolder\n");

  #We check if the config file is already present in the given folder. If this happens, the user can decide the overwriting
  if ( -e ($sessionFolder.'/'.$ARGV[0]) ){
      my $sure= get_input("A config file $ARGV[0] is already present in the folder $sessionFolder. Do you want to overwrite it?(y or n) [Suggested: y]",'^[YynN]$');
          annoPrint($sure."\n",1);
          if ( ( $sure eq "y") or ( $sure eq "Y")){
            (system("cp $ARGV[0] $sessionFolder") ) == 0 
              or annoDie("ERROR: Unable to copy $ARGV[0] in $sessionFolder. Please check permissions...\n");
          }
          else { annoPrint(	"Ok, Annocript will use the one present there...\n");}
  }else{
    (system("cp $ARGV[0] $sessionFolder") ) == 0 
        or annoDie("ERROR: Unable to copy $ARGV[0] in $sessionFolder. Please check permissions...\n");
  }
}



=head2 annoPrint

 Title   : annoPrint
 Usage   : annoPrint( - string -> the sentence that have to be annoPrint ();

 Function: will print (the string in input always in the log file and on the STDOUT
 
 Returns : nothing

=cut
sub annoPrint{
  my $string = shift;    
  my $onlyLog = shift;
  
  if ( defined $onlyLog){
    print $logFHandle $string;
  }else{
    my $STDOUT = *STDOUT;
    #Prints on both OUT
    for ($logFHandle, $STDOUT) { print $_ $string; }
  }
}


=head2 annoDie

 Title   : annoDie
 Usage   : annoDie( - string -> the sentence that have to be print

 Function: will print the string and die
 
 Returns : nothing

=cut
sub annoDie{
  my $string = shift; 
       
  my $STDOUT = *STDOUT;
  
  #Print on both the log and the STDOUT the error
  for ($logFHandle, $STDOUT) { print $_ $string; }
  
  if ( -d $logFolder){
    #Moving the log file to the log folder before to exit Annocript
    move($workingFolder.'/'.$logFile, $logFolder."/".$logFile) or die "Cannot move $workingFolder/$logFile to $logFolder/$logFile\n";
  }else{
    print "Annocript will not move the log file: $workingFolder/$logFile. All the log files in $workingFolder will be removed".
    " at next successful Annocript run.\n";
  }
  #Die and print error number
  die;# $!;
}


=head2 nice_printing

 Title   : nice_printing
 Usage   : nice_printing( - sentence -> the sentence that have to be annoPrint ();

 Function: this subroutine prints in a box of asterisks the sentence in input
 
 Returns : nothing

=cut
sub nice_printing{
    my $sentence = shift;
    
    my $length = length($sentence);
    my $i;
    for ($i=0; $i<($length);$i++){
      annoPrint ("#");
    }
    annoPrint ("####\n# $sentence #\n####");
    for ( $i=0; $i<$length; $i++){
		annoPrint ("#");
    }
    annoPrint ("\n\n");
}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configUser -> file with the user configuration
                              - configAnnocript -> file with the basic parameters of Annocript
                              - variablesFile -> the path to a file with all variables written
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.
 
 Returns : nothing

=cut
sub checkConfigVariablesOLD {
  my $configUser = shift;
  my $configAnnocript = shift;
  my $variablesFile = shift;
  
  my $hashCheck;        
  if (! open(VARF,"<$variablesFile")){ annoDie ("ERROR: Failure opening '$variablesFile'. Your Annocript version is corrupted - $!");}
  if (! open(CUSER,"<$configUser")){ annoDie ("ERROR: Cannot find '$configUser' - You can pick it from Annocript folder and  put in $workingFolder.");}
  if (! open(CANN,"<$configAnnocript")){ annoDie ("ERROR: Failure opening '$configAnnocript'. Your Annocript version is corrupted - $!");}
  
  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CUSER>){ 
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashCheck->{$1} = "OK";
    }
  }	
  #Stores the variables in the config annocript file inside the hash
  $start = 0;
  while (my $line = <CANN>){ 
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
      $hashCheck->{$1} = "OK"; 
      #annoPrint ($line."\n"); #DEBUGCODE
    }
  }	
  close(CUSER);
  close(CANN);
    
  my @allVars = ();  
  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $line=0;
  #For each line of the variables file
  while (my $line = <VARF>){ 
    
    #get the variables in the line
    my @variables = split (/;/,$line);
    
    #put the variables inside an array
    push (@allVars, @variables);
    
    $line++;
    #For each of the variables in the variables file
    foreach my $var (@variables){
      #annoPrint ("Variable: $var - - value: ".$hashCheck->{$var}."\n");#DEBUGCODE
      $var =~ s/\n//;
      if( !(defined($hashCheck->{$var})) ){
          if ($line == 1){
            annoDie( "ERROR: in $configAnnocript variable $var is missing. Please check the file. Closing Annocript...\n "); 
            $errors=1;
          }
          if ($line == 2){
            annoDie( "ERROR: in $ARGV[0] - variable $var is missing. Please check the file. Closing Annocript...\n "); 
            $errors=1;
          }
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
    }
    
  }
  
  #print_array(\@allVars);	
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashCheck){
     # print "Search $key in array...\n";#DEBUGCODE
      if (!(grep {/$key/} @allVars )){
          annoDie("ERROR: Variable $key is in the config files and not in $variablesFile file. This is completely wrong. Re-install Annocript.\n ");
        }
    }
  
  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configUser -> file with the user configuration
                              - configAnnocript -> file with the basic parameters of Annocript
                              - variablesFile -> the path to a file with all variables written
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.
 
 Returns : nothing

=cut
sub checkConfigVariables {
  my $configUser = shift;
  my $configAnnocript = shift;
  my $variablesFile = shift;
  
  my $hashUserCheck;
  my $hashAnnCheck;
          
  if (! open(VARF,"<$variablesFile")){ annoDie ("ERROR: Failure opening '$variablesFile'. Your Annocript version is corrupted - $!");}
  if (! open(CUSER,"<$configUser")){ annoDie ("ERROR: Cannot find '$configUser' - You can pick it from Annocript folder and  put in $workingFolder.");}
  if (! open(CANN,"<$configAnnocript")){ annoDie ("ERROR: Failure opening '$configAnnocript'. Your Annocript version is corrupted - $!");}
  
  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CUSER>){ 
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashUserCheck->{$1} = "OK";
    }
  }	
  #Stores the variables in the config annocript file inside the hash
  $start = 0;
  while (my $line = <CANN>){ 
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
      $hashAnnCheck->{$1} = "OK"; 
      #annoPrint ($line."\n"); #DEBUGCODE
    }
  }	
  close(CUSER);
  close(CANN);
    
  my @userConfVars = ();  
  my @annConfVars = ();
  
  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $lines=0;
  #For each line of the variables file
  while (my $line = <VARF>){ 
    $line =~ s/\n//;#Remove \n in the end of the line
    
    #get the variables in the line
    my @variables = split (/;/,$line);
    
    $lines++;
    #For each of the variables in the variables file
    foreach my $var (@variables){
      
      #$var =~ s/\n//;#Remove
      
      if ($lines == 1){
        #annoPrint ("Annocript Variable: $var - - value: ".$hashAnnCheck->{$var}."\n");#DEBUGCODE
        #put the variable inside an array for ann config
        push (@annConfVars, $var);
        if( !(defined($hashAnnCheck->{$var})) ){
          
          annoDie( "ERROR: in $configAnnocript variable $var is missing. Please check the file. Closing Annocript...\n "); 
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      }
      if ($lines == 2){
        #annoPrint ("User Variable: $var - - value: ".$hashUserCheck->{$var}."\n");#DEBUGCODE
        #put the variable inside an array for user config
        push (@userConfVars, $var);
        if( !(defined($hashUserCheck->{$var})) ){
          annoDie( "ERROR: in $ARGV[0] - variable $var is missing. Please check the file. Closing Annocript...\n "); 
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      } 
    }
  }
  
  #print_array(\@allVars);	
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashUserCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @userConfVars )){
       annoDie("ERROR: Variable $key is in the user config files and not in $variablesFile file. This is completely wrong. Annocript will not work...\n ");
     }
  }
  
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashAnnCheck){
    # print "Search $key in array...\n";#DEBUGCODE
    if (!(grep {/$key/} @annConfVars )){
      annoDie("ERROR: Variable $key is in the annocript config files and not in $variablesFile file. This is completely wrong. Annocript will not work...\n ");
    }
  }

  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}

=head2 checkVariable

 Title   : checkVariable
 Usage   : checkVariable(  - var -> value of the variable to check
                           - name -> name of the variable
                           - sentence -> something to write to the user in case of error)

 Function:  this subroutine checks if a variable has the YES or NO value. Dies otherwise.
 
 Returns : nothing

=cut
sub checkVariable {
  my $var = shift;
  my $name = shift;
  my $yesSentence = shift;
  my $noSentence = shift;
  
  my $no = '';
  if (defined $noSentence){
      $no = $noSentence;
  }
  if ( ($var ne 'YES') and ($var ne 'NO')){
    annoDie("ERROR: Check variable $name in the config file. $var is a wrong value!\n");
  }elsif ($var eq 'YES'){
    annoPrint ($yesSentence);
  }elsif ($var eq 'NO'){
    annoPrint ($no);
  }
}


  
=head2 createBasicFolders

 Title   : createBasicFolders
 Usage   : createBasicFolders(  );

 Function:  this subroutine creates the basic folders for Annocript. If the folders.txt file is not present 
            means that Annocript has not been installed in the folder used.
 
 Returns : nothing

=cut
sub createBasicFolders{

  #If the file with folders references is installed Annocript can start
  open (FOLD, $foldersFile) or annoDie(" ERROR: Annocript has not been installed here. Please execute install.pl in the program folder.\n");
  my $line = <FOLD>;
  
  #Extracts folders paths
  my @folders = split(" ", $line);
  
  #The working folder is the first path
  $workingFolder = $folders[0];
  
  #Cancel the final slash..if it is there it is removed
  if($workingFolder =~ /\/$/){
    chop($workingFolder);
  }
  
  $mainDataFolder = $workingFolder."/".$mainDataFolder;
  $jobsFolder = $workingFolder."/".$jobsFolder;
  $sessionFolder = $jobsFolder;
  $annocriptFolder = $folders[1];
  #$workingFolder = getcwd;
  
  $configHash->{'dataFolder'} = $dataFolder;
  $configHash->{'gffFolder'} = $gffFolder;
  $configHash->{'outFolder'} = $outFolder;
  $configHash->{'statsFolder'} = $statsFolder;
  $configHash->{'logFolder'} = $logFolder;
  close(FOLD);
  
  #Creates the jobs folder if it does not exists
  unless(-d $jobsFolder){
    annoPrint ($jobsFolder." doesn't exists. Creating folder...\n");
    mkdir $jobsFolder or annoDie("ERROR: can't create folder $jobsFolder\n ");
  }
   #Creates the main 'data' folder if it does not exists  
  unless(-d $mainDataFolder){
    annoPrint ($mainDataFolder." doesn't exists. Creating folder...\n");
    mkdir $mainDataFolder or annoDie( "ERROR: can't create folder $mainDataFolder\n ");
  }
  
  #annoPrint ("Annocript folder: $annocriptFolder - Users: $jobsFolder - \n");
  
}




=head2 showOptions

 Title   : showOptions
 Usage   : showOptions(  );

 Function:  Shows some options the user can choose among.
 
 Returns : 1 if the user choose a new session else 0

=cut	
sub showOptions{
  
  my $newSession = 0;
  
  my $myFolder;
  annoPrint ("1. New session - 2. Use Previous - 3. Remove Session - 4. Do some utils - 5. Exit\n");
  my $session = get_input("Your choice: ","^[12345]\$");
   annoPrint ($session."\n",1);#Print on the log file
   
  #A  construct will execute different code depending from the choose
    if ($session eq "1") { $newSession = 1; $myFolder = changeSession();}
	  if ($session eq "2") { 
	      if ( !(is_folder_empty($jobsFolder)) ){
	        my $workDir = getcwd;
	        chdir $jobsFolder;
	        (system("ls -d */")) == 0 or annoDie ("Unable to list directory $jobsFolder!\n"); 
	        annoPrint ("Choose your folder: ");
	        $myFolder = <STDIN>;
	        chomp $myFolder;
	        $myFolder = $1 if($myFolder=~/(.*)\/$/);#Cut last letter if it is /
	        #print $myFolder;
	        my $validSession = 0;
	        if (valid_session($myFolder) == 1){
							$validSession = 1;
					}
	        while ( !(-d $myFolder)){
	          (system("ls -d */")) == 0 or annoDie ("Unable to list directory $jobsFolder!\n"); 
	          annoPrint ("Choose your folder: ");
	          $myFolder = <STDIN>;
	          chomp $myFolder;
	          $myFolder = $1 if($myFolder=~/(.*)\/$/);#Cut last letter if it is /
	          if (valid_session($myFolder) == 1){
							$validSession = 1;
						}
	        }
	        $sessionFolder =  $sessionFolder."/".$myFolder;
	        chdir $workDir;      
	      }else{
	        annoPrint("There are not existing sessions. Please restart Annocript and create a new one!\n");
	        exit 1;
	      }
    }
    if ($session eq "3") {
      my $workDir = getcwd;
      chdir $jobsFolder;
      (system("ls")) == 0 
        or annoDie ("Error: Unable to list directory $jobsFolder. Check permissions!\n"); 
      annoPrint ("Choose the session to remove: ");
      $myFolder = <STDIN>;
      chomp $myFolder;
			while(!(-d $myFolder)){
        (system("ls")) == 0 
        or annoDie ("Error: Unable to list directory $jobsFolder. Check permissions!\n"); 
				annoPrint ("Choose your folder: ");
				$myFolder = <STDIN>;
				chomp $myFolder;
      }
      my $sessCancel = get_input("Really you want to remove $myFolder ?(y or n)",'^[YynN]$');
      #Code when the user wants to overwrite
      annoPrint($sessCancel."\n",1);#Print on the log file
      if ( ($sessCancel eq "y") or ($sessCancel eq "Y") ){
        my $sure= get_input("Are you sure?(y or n)",'^[YynN]$');
        annoPrint($sure."\n",1);#Print on the log file
          if ( ( $sure eq "y") or ( $sure eq "Y")){
            remove_tree($myFolder);
            print "$myFolder correctly removed!\n";
            exit 1;
          }
          else { print	"$myFolder will not be removed!\n"; exit 1;}
      }else{ print "$myFolder will not be removed!\n"; exit 1;} 
      #chdir $workDir;
      
     }
    if ($session eq "4"){
      my $configAnnPath = $annocriptFolder."/".$configAnnocript;
      my $script = $annocriptFolder."/".$utilsScript;
      #(system("perl $script $configAnnPath")) == 0 
        #or annoDie("Error: Unable to call $script!\n");
      execute_utils($configAnnPath);
      exit 1;
      } 
   if ($session eq "5")  {annoPrint ("Annocript was interrupted by user... bye!\n");exit 1;}

  annoPrint ($myFolder."\n",1);
  
 return $newSession;
}

=head2 valid_session

 Title   : valid_session
 Usage   : valid_session( $session -> the name of the folder );

 Function:  When user wants to open an existing session Annocript will controll if this does not come
            from a wrong execution: some folders must be inside
 
 Returns : 1 if the session is valid, else -1

=cut
sub valid_session{
  my $session = shift;
  
  my $missFolders = "";
  my $retVal = -1;
  
  if (-d $session){		
		#Folders to be alive in the session folder if it is consistent
		if ( !(-d $session."/".$dataFolder) ){
			$missFolders.= " ".$dataFolder;
		}
		
		if ( !(-d $session."/".$logFolder) ){
			$missFolders.= " ".$logFolder; 
		}
		
		if ( !(-d $session."/".$outFolder) ){
			$missFolders.= " ".$outFolder;
		}
		
		if ( !(-d $session."/".$statsFolder) ){
			$missFolders.= " ".$statsFolder;
		}
		
		if ( !(-d $session."/".$gffFolder) ){
			$missFolders.= " ".$gffFolder;
		}
		
		if ( $missFolders eq ""){
				$retVal = 1;
		}else{
			my $sessCancel = get_input("Session $session is not valid. Do you want to remove it?(y or n)",'^[YynN]$');
			
			#Code when the user wants to overwrite
			annoPrint($sessCancel."\n",1);#Print on the log file
			if ( ($sessCancel eq "y") or ($sessCancel eq "Y") ){
				my $sure= get_input("Are you sure?(y or n)",'^[YynN]$');
				annoPrint($sure."\n",1);#Print on the log file
					if ( ( $sure eq "y") or ( $sure eq "Y")){
						remove_tree($session);
						annoPrint("Done! Annocript will exit..\n");
						exit 1;
					}
				 # else { annoPrint(	"Ok, choose another name. But don't let Annocript waste its time...\n");}
			}
			#else { annoPrint("Ok, choose another name. But don't let Annocript waste its time...\n");}
			
		#print "Session $session is not valid. Please remove it.\n";
		#
		print "Error: folders $missFolders are missing\n";
		}
 }else{print "Folder $session does not exist.\n";}
  return $retVal;
}



=head2 changeSession

 Title   : changeSession
 Usage   : changeSession( );

 Function:  When user wants to start another session it asks about the name and check if that name is already been taken. At the end
			of the function the folder is created.
 
 Returns : nothing

=cut
sub changeSession{
  annoPrint ("Choose a name for your session: "); 
  my $session = <STDIN>;
  chomp $session;
	
  
  #Check that the use does not insert strange chars in the name
  while ( $session !~ /^[A-za-z0-9\_\-]*$/i ){
		annoPrint( "The name $session is not permitted. \nChoose another using [A-za-z0-9\_\-]:");
    $session = <STDIN>;
    chomp $session;
	}   
   
  
  #Controls if the session folder already exists and does something
  while (-d $sessionFolder."/".$session){
    my $ok = 0;
    my $sesNaOverWrite = get_input("This name has already been taken, overwrite it?(y or n)",'^[YynN]$');
    
    #Code when the user wants to overwrite
    annoPrint($sesNaOverWrite."\n",1);#Print on the log file
    if ( ($sesNaOverWrite eq "y") or ($sesNaOverWrite eq "Y") ){
      my $sure= get_input("Are you sure?(y or n)",'^[YynN]$');
      annoPrint($sure."\n",1);#Print on the log file
        if ( ( $sure eq "y") or ( $sure eq "Y")){
          remove_tree($sessionFolder."/".$session); 
          $ok = 1; 
        }else { annoPrint(	"Ok, choose another name. But don't let Annocript waste its time...\n");}
    }else { annoPrint("Ok, choose another name. But don't let Annocript waste its time...\n");}	
		
    #If the user wants to choose a new name
    if($ok == 0){
      annoPrint ("Choose a name for your session: ");
      $session = <STDIN>;
      chomp $session;
    }
  } 
  #Assigns the path to a global variable
  $sessionFolder = $sessionFolder."/".$session;
  
  #Creates the folder
  mkdir $sessionFolder or annoDie("ERROR: Unable to create the directory... please check your permissions.Exiting..\n");
  annoPrint ($sessionFolder." created!\n");
  
  return $session;
}
	 

=head2 get_input

 Title   : get_input
 Usage   : get_input( - sentence: a sentence that will printed in input to ask something to the user;
					 - regex: what the answer of the user have to respect
                               );

 Function: Takes in input a sentence and a regex. Asks to the user the sentence and controls its input with regex
 
 Returns : input given by the user

=cut
sub get_input{
  my $sentence = shift;
  my $regex = shift;
  my $input='';
		
  while (!($input =~ /$regex/) or ($input eq "")){  
    annoPrint ($sentence." ");
    $input = <STDIN>;
    chomp $input;
    annoPrint ("\n");
  }	
  return $input;
}



=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
 Returns : nothing

=cut
sub configFile2Hash{  
  my $configFilePath=shift;
  
  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or annoDie("ERROR: The file $configFilePath doesn't exists. Annocript will exit..\n");
  while (my $line = <configFile>){ 
    if ($line =~ /#########/){
      $start = 1;
    }
    if( ($line =~ /(\S+)\s*=\s*(\S+)/) and ($start == 1) and !($line =~ /#/) ){
      $configHash->{$1} = $2; 
      #annoPrint ("$1 = $2\n") ;#DEBUGCODE     
    }
  }	 
  close(configFile);	
	#annoPrint (Dumper\$configHash); #DEBUGCODE
}


=head2 completeConfiguration

 Title   : completeConfiguration
 Usage   : completeConfiguration( - configFilePath = path of the config file
                               );

 Function:  this subroutine completes the configuration of annocript:
              - checks the variables with YES, NO values
              - creates needed folders
              - check variables related to the three modules (DB_CREATOR, PROGRAMS_EXEC, GFF3_AND_OUTPUT)
              - build some important variables in the configuration hash
 Returns : nothing

=cut
sub completeConfiguration{
   
  #Here I put a control of the variables in the configuration file that can be YES or NO
  #annoPrint ("\n \n");
  checkVariable($configHash->{'doDbCreation'},'doDbCreation', "DB file will be downloaded from UNIPROT website and DB built.\n");
  checkVariable($configHash->{'doExecutePrograms'},'doExecutePrograms', "Some programs will be run. In the following are listed\n");
  
  if ($configHash->{'doExecutePrograms'} eq 'YES'){
    checkVariable($configHash->{'doBlastx'},'doBlastx',"A BLASTX of your sequences will be done.\n");
    checkVariable($configHash->{'doRpstblastn'},'doRpstblastn',"A RPS BLAST of your sequences will be done.\n");
    checkVariable($configHash->{'doBlastn'},'doBlastn',"A BLASTN of your sequences will be done.\n");
    checkVariable($configHash->{'doPortrait'},'doPortrait',"A PORTRAIT of your sequences will be done.\n");
    checkVariable($configHash->{'doDna2Pep'},'doDna2Pep',"DNA2PEP on your sequences will be done.\n");
  }
  checkVariable($configHash->{'doBuildOutput'},'doBuildOutput', "The final table will be built\n");
  checkVariable($configHash->{'extractStatistics'},'extractStatistics',"A web page and file with statistics will be created\n");
  
  checkVariable($configHash->{'useGFFDB'},'useGFFDB',"GFF3 files will be printed\n");
  
  
  #The configuration is shown and the user should accept or not
  my $agreeWithConfig = get_input("Do you agree with it?(y or n)",'^[YynN]$');
  annoPrint($agreeWithConfig."\n",1 );
    if ( ($agreeWithConfig eq "y") or ($agreeWithConfig eq "Y") ){
      annoPrint ("Ok!.. let's go!\n");
    }else{
      annoPrint("Ok! Let's change the configuration file! Annocript will close...\n");
      exit 1;
      }
  
  
	#Storing the session folder in the hash  
  $configHash->{'sessionFolder'} = $sessionFolder;
	
  #Checks everything regarding the query file. If session is new, the file is put in the session folder
  check_query_file();
  
  $dataFolder = $sessionFolder."/".$configHash-> {'dataFolder'};
  #Check if directory DATA exists, otherwise it creates it
  unless(-d $dataFolder){
    annoPrint ("Creating folder $dataFolder...\n");
    mkdir $dataFolder or annoDie("ERROR: can't create folder $dataFolder. Check permissions. \n");
  }
  
  $logFolder = $sessionFolder."/".$configHash-> {'logFolder'};
  #Check if directory LOG exists, otherwise it creates it
  unless(-d $logFolder){
    annoPrint ("Creating folder $logFolder...\n");
    mkdir $logFolder or annoDie("ERROR: can't create folder $logFolder. Check permissions. \n");
  }
  
  $outFolder = $sessionFolder."/".$configHash->{'outFolder'};
  #Check if directory OUTPUT exists, otherwise it creates it
  unless(-d $outFolder){
    annoPrint ("Creating folder $outFolder....\n");
    mkdir $outFolder or annoDie("ERROR: can't create folder $outFolder. Check permissions. \n");
  }
  
  $statsFolder = $sessionFolder."/".$configHash->{'statsFolder'};#Statistics folder
  #Creation of a folder for the statistics
  unless(-d $statsFolder){
    annoPrint ("Creating folder $statsFolder...\n");
    mkdir $statsFolder or annoDie("ERROR: can't create folder $statsFolder. Check permissions. \n");
  }#else{annoPrint ("Folder $statsFolder already exists. Continuing...\n ";}
 
  $gffFolder = $sessionFolder."/".$configHash->{'gffFolder'};
  #Check if directory for GFF files exists, otherwise it creates it
  unless(-d $gffFolder){
    annoPrint ("Creating folder ".$gffFolder."...\n");
    mkdir $gffFolder or annoDie( "ERROR : can't create folder $gffFolder. Check permissions. \n");
  } 
 
  #Check if directory for BLASTX output exists, otherwise it creates it
  unless(-d $outFolder."/blastx"){
    annoPrint ("Creating folder ".$outFolder."/blastx...\n");
    mkdir $outFolder."/blastx" or annoDie( "ERROR : can't create folder $outFolder /blastx. Check permissions. \n");
  } 
  #Check if directory for RPSTBLASTN output exists, otherwise it creates it
  unless(-d $outFolder."/rpstblastn"){
    annoPrint ("Creating folder ".$outFolder."/rpstblastn...\n");
    mkdir $outFolder."/rpstblastn" or annoDie( "ERROR : can't create folder $outFolder /rpstblastn . Check permissions. \n");
  }  
  #Check if directory for BLASTN output exists, otherwise it creates it
  unless(-d $outFolder."/blastn"){
    annoPrint ("Creating folder ".$outFolder."/blastn...\n");
    mkdir $outFolder."/blastn" or annoDie( "ERROR : can't create folder $outFolder /blastn. Check permissions. \n");
  }
   
      
  #annoPrint (Dumper\$configHash;    #DEBUGCODE

  
  annoPrint ("\n>>> Database settings (it will need a while):\n");
	checkDB_CREATION($configHash-> {'DBCreationFolder'});
	annoPrint ("\n>>> Programs settings:\n");
	checkPROGRAMS_EXEC($configHash->{'ProgExecFolder'});
	annoPrint ("\n>>> Setting final statistics and outputs\n");
	checkGFF3_AND_OUTPUT($configHash-> {'GFFAndOutputFolder'});
	annoPrint ("..OK! WELL DONE!\n");
	
	my $dbDataFolder = $configHash->{'dbDataFolder'};
	
  #Adding some variables for blast out types. This is done 
  #because if a day we want to add more types, nothing is to change in gff3_manager.pm
  $configHash->{'blastxSPOutType'} = $configHash->{'blastOutType'};
  $configHash->{'blastxTROutType'} = $configHash->{'blastOutType'};
  $configHash->{'blastxUnirefOutType'} = $configHash->{'blastOutType'};
  $configHash->{'rpstblastnOutType'} = $configHash->{'blastOutType'};
  $configHash->{'blastnOutType'} = $configHash->{'blastOutType'};
  
  #Writing the complete path (with datafolder) of the databases    
  $configHash->{'rfamDB'} = "$dbDataFolder/".$configHash->{'rfamDB'};
	
  if ( $configHash->{'dbInUse'} eq "uniprotkb" ){
    $configHash->{'swissProtDB'} = "$dbDataFolder/".$configHash->{'swissProtDB'};
    $configHash->{'tremblDB'}    = "$dbDataFolder/".$configHash->{'tremblDB'};
    $configHash->{'uniprotKBDB'} = "$dbDataFolder/".$configHash->{'uniprotKBDB'};
  }elsif (  $configHash->{'dbInUse'} eq "uniref" ){
    $configHash->{'swissProtDB'} = "$dbDataFolder/".$configHash->{'swissProtDB'};
    $configHash->{'unirefDB'} = "$dbDataFolder/".$configHash->{'unirefDB'};
  }
 
  #Now give to the output the names with the input fasta file name
  my $preamble = $configHash->{'name'}.'_'.$configHash->{'dbInUse'}.'_'.$configHash->{'dbVersion'}.'_';
  $configHash->{'outFileName'} = $preamble.$configHash->{'outFileName'};
  $configHash->{'outFiltered'} = $preamble.$configHash->{'outFiltered'};
  $configHash->{'ORFFastaFileName'} = $preamble.$configHash->{'ORFFastaFileName'};
  $configHash->{'NCOutFile'} = $preamble.$configHash->{'NCOutFile'};
  $configHash->{'codingOutFile'} = $preamble.$configHash->{'codingOutFile'};

  #Setting paths for the stats image folder and the html pages
  $configHash->{'imgDir'} = $preamble."statspage";  
  $configHash->{'imgPath'} = $statsFolder."/".$configHash->{'imgDir'};
  
  #The version of Annocript
  $configHash->{'AnnocriptVer'} = $version;
  $configHash->{'dbTablesVersion'} = $dbTablesVersion;
  $configHash->{'usefulFolder'} = $usefulFolder;
  $configHash->{'faSomeRecords'} = $annocriptFolder."/$usefulFolder/".$configHash->{'faSomeRecords'}.
  check_dependencies();
}



=head2 detect_gff_db

 Title   : detect_gff_db
 Usage   : detect_gff_db( -db => 'name of the database,
						-platform, -host, -port, -mySqlUser, -mySqlPass => 'Parameters for DB using and opening'
                               );

 Function:  check if a GFF database exists. If it does, then it will ask if the user wants to exit Annocript
            otherwise the database will be removed.

 Returns : nothing

=cut
sub detect_gff_db{
  
  
  my $dbToSearch = shift;
	my $platform = shift;
	my $host = shift;
	my $port = shift;
	my $mySqlUser = shift;
	my $mySqlPass = shift;
  
  # PERL DBI CONNECT AND CREATE TABLE
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } );
  my $databases = $dbh->selectcol_arrayref("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '$dbToSearch'");
  
  my $numKeys = scalar(@$databases);
  
  #The db exists
  if( $numKeys >0 ){
    my $removeDb = get_input("Existent db $dbToSearch. If you chose to obtain the final output, it'll be removed. Continue?(y or n)[Suggested: Y]\n",'^[YynN]$');
     annoPrint($removeDb."\n",1);#Print on the log file
    if ( ($removeDb eq "y") or ($removeDb eq "Y") ){
      my $sure= get_input("Are you sure?(y or n)",'^[YynN]$');
      annoPrint($sure."\n",1);
      if ( ($sure eq "y") or ($sure eq "Y") ){
				$configHash->{'destroy_gffdb'} = "YES";  
        annoPrint ("Annocript will continue....\n"); 
      }
      else { 
		    $configHash->{'destroy_gffdb'} = "NO"; 
        annoPrint( "Ok, choose another name in the config file and don't let Annocript waste its time...\n");
        exit 1;
      }	
    }
    #The db exists and the user doesn't wants to cancel it			
    else { 
      annoDie( "Ok, choose another name in the config file and don't let Annocript waste its time...\n");
      exit 1;
    }	
    $dbh->disconnect();
  }else{$configHash->{'destroy_gffdb'} = "NO"; }
}






=head2 uniprot_db_selection

 Title   : uniprot_db_selection
 Usage   : uniprot_db_selection( -db => 'name of the database,
						-platform, -host, -port, -mySqlUser, -mySqlPass => 'Parameters for DB using and opening'
                               );

 Function:  check if a database exists. It access to information_schema database and uses the
            show database SQL command to create a list where one can find the specific database.

 Returns : nothing

=cut
sub uniprot_db_selection{
  my $dbToSearch = shift;
  my $platform = shift;
  my $host = shift;
  my $port = shift;
  my $mySqlUser = shift;
  my $mySqlPass = shift;
  my $doDbCreation = shift;
  my $doBuildOutput =shift;
	
  my $dsn = "dbi:$platform:information_schema:$host:$port";
              
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1} );
  my $databases = $dbh->selectcol_arrayref("show databases like '".$dbToSearch."'");
  
  #Get the number of results of the LIKE SQL operator
  my $numKeys = @$databases;
 
  #The db has to be created
  if( $doDbCreation eq "YES"){
    #The db exists
    if( $numKeys >0 ){
      #Now that the user wanted to overwrite an existing database, 
      #he can do it only if the database has the same coding version of the one given in this script
      #DON'T FORGET TO UPDATE THE VARIABLE (dbTablesVersion) EACH TIME THE DB STRUCTURE WILL CHANGE!!!
      my $dsn2 = "dbi:$platform:$dbToSearch:$host:$port";
      check_db_tables_version($dsn2,$mySqlUser,$mySqlPass) 
		or annoDie( "In this version of Annocript ($version) the database will have different structure.". 
		"You cannot update a previous version. Please create a completely new database.\n");
      
      my $removeDb = get_input("$dbToSearch already exists. It will be checked and completed if necessary. Continue?(y or n)\n",'^[YynN]$');
      annoPrint($removeDb."\n",1);
      if ( ($removeDb eq "y") or ($removeDb eq "Y") ){
        my $sure= get_input("Are you sure?(y or n)",'^[YynN]$');
        annoPrint($sure."\n",1);
        if ( ($sure eq "y") or ($sure eq "Y") ){
          #$dbh->do("DROP DATABASE IF EXISTS $dbToSearch") or die "ERROR: unable to cancel database $dbToSearch!"; 
          #IN THIS CASE ANNOCRIPT WILL DO NOTHING. Later in the db_creator module it will remove the database.
          #You could also remove the database here but if a programmer wants to use Annocript to write only one table in the database
          #he can't. To do this now simply comment the createDBandTables subroutine execution.
          $configHash->{'dbExistent'} = 'YES';
        }
        else { die "Ok, please pay attention with the database! Annocript is exiting...\n";}	
      }
      #The db exists and the user doesn't wants to cancel it			
      #else { print	"Database $db will be removed if you choosed in the config file the CreateDB=YES options...\n";}	
      else { die "Ok, please pay attention with the database! Annocript is exiting...\n";}
      $dbh->disconnect();
    }else{$configHash->{'dbExistent'} = 'NO';}
  }
  #Do not create database for proteins but it already exists and we want to use to build output
  elsif(( $doDbCreation eq "NO") and ( $doBuildOutput eq "YES") ){
    #If it doesn't exists than user has to reconfigure Annocript
    if( $numKeys == 0 ){
      annoDie( "ERROR: Database $dbToSearch you choosed in the config file doesn't exists. Please check it. Annocript will close...\n");
    }
    #if exists Annocript will continue
    else{
      $dbh->disconnect();
      $configHash->{'dbExistent'} = 'YES';
    }
  }else{
    #We don't have to cancel the database. Simply disconnect it if it exists
    if( $numKeys >0 ){
      $dbh->disconnect();
      $configHash->{'dbExistent'} = 'YES';
    }
  }
	return 1;
}

=head2 check_db_tables_version

 Title   : check_db_tables_version
 Usage   : check_db_tables_version( - database = the db to check in
                               );

 Function:   says if the version given in input is corresponding to the one of the database going to update
 Returns : 1 if ok, else 0

=cut
sub check_db_tables_version {
  my $dsn = shift;
  my $mySqlUser = shift;
	my $mySqlPass = shift;
  
  my $retVal = 0;
  
  my $dbh = DBI->connect($dsn,$mySqlUser,$mySqlPass, { RaiseError => 1} );
  
  my $select_handle = $dbh->prepare("SELECT * FROM infoTable WHERE dbTablesVersion='".$dbTablesVersion."'");
  die "ERROR [$?]: can't prepare queries; aborting: ?\n" unless defined $select_handle;
  $select_handle->execute() or die $select_handle->errstr;
  my @res = $select_handle->fetchrow_array;
  if (scalar(@res) != 0){
     $retVal = 1;
  }
  return $retVal;
}


=head2 checkDB_CREATION

 Title   : checkDB_CREATION
 Usage   : checkDB_CREATION( - configHash = the piece of the config hash for DB_CREATION
							 - local folder = the folder where actually is the DB_CREATION script
                               );

 Function:   Checks the config hash parameters regardin the database creation and user. How they are written, if they are present and so on..
 Returns : nothing

=cut
sub checkDB_CREATION {
  my $localFolder = shift;

  my $wrongLinks = 0;  
 
  ##########################CHECKING THE LINKS   
  if ($configHash->{'doDbCreation'} eq 'YES'){
    
    #The first check is on the permissions. If the user does not have it cannot create database or perform analysis
    if (!check_DB_user_permissions($configHash->{'platform'},$configHash->{'host'},
              $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'})) {
        annoDie("ERROR: You do not have enough privileges to create the DB. Please ask your adiministrator. Annocript will close...\n");
    }
      # annoPrint (Dumper\$configHash;
      
      #A check on email address 
      if ( !($configHash->{'uniprotWebPass'} =~ /@/) or !($configHash->{'uniprotWebPass'} =~ /\.[a-zA-Z]+?/) ){
        annoDie( "ERROR : ".$configHash->{'uniprotWebPass'}." is not a correct email address. Please check it in config file \n");
      }
      
      #Triple controls on the links. We need it three times because it can fail.
      annoPrint ("Checking the links for the databases to download...(it can take a while)...");
			
      #This IFs contain the check of the links to the databases. Two types of check are present for links to the Uniprot data.
      # LWP::Simple::Head and NET::Ftp
      #print "Swiss-prot DB link...";
        if (!head($configHash->{'swissprotDBLink'}) ) {
          if (testFTP($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              my $new_link = check_FTP_diff_sources($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
              if ( $new_link ne ''){
                  $configHash->{'swissprotDBLink'} = $new_link;
              }else{
                annoPrint ("Please check if [".$configHash->{'swissprotDBLink'}."] is a correct URL. Annocript will continue...\n"); #
                $wrongLinks++;
              }
            }
          }
        }
       #print "Trembl DB link...";
        if (!head($configHash->{'tremblDBLink'}) ) {
          if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              my $new_link = check_FTP_diff_sources($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                if ( $new_link ne ''){
                    $configHash->{'tremblDBLink'} = $new_link;
                }else{
                  annoPrint ("Please check if [".$configHash->{'tremblDBLink'}."] is a correct URL. Annocript will continue...\n"); 
                  $wrongLinks++;
                }
              }
            }
          }
        }
        #print "Uniprot DB version link...";
        if (!head($configHash->{'uniprotVerLink'}) ) {
          if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                my $new_link = check_FTP_diff_sources($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                if ( $new_link ne ''){
                    $configHash->{'uniprotVerLink'} = $new_link;
                }else{
                  annoPrint ("Please check if [".$configHash->{'uniprotVerLink'}."] is a correct URL. Annocript will continue...\n"); #
                  $wrongLinks++;
                }
              }
            }
          }
        }  
        #print "UniRef DB link...";
        if (!head($configHash->{'unirefDBLink'}) ) {
          if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                my $new_link = check_FTP_diff_sources($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                if ( $new_link ne ''){
                    $configHash->{'unirefDBLink'} = $new_link;
                }else{
                  annoPrint ("Please check if [".$configHash->{'unirefDBLink'}."] is a correct URL. Annocript will continue...\n"); #
                  $wrongLinks++;
                }
              }
            }
          }
        }
        #print "Uniref DB version link...";
        if (!head($configHash->{'unirefVerLink'}) ) {
          if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                my $new_link = check_FTP_diff_sources($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                if ( $new_link ne ''){
                    $configHash->{'unirefVerLink'} = $new_link;
                }else{
                  annoPrint ("Please check if [".$configHash->{'unirefVerLink'}."] is a correct URL. Annocript will continue...\n"); #
                  $wrongLinks++;
                }
              }
            }
          }
        }
        #print "Uniprot idmapping link...";
        if (!head($configHash->{'GODBLink'}) ) {
          if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                my $new_link = check_FTP_diff_sources($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                if ( $new_link ne ''){
                    $configHash->{'GODBLink'} = $new_link;
                }else{
                  annoPrint ("Please check if [".$configHash->{'GODBLink'}."] is a correct URL. Annocript will continue...\n"); 
                  $wrongLinks++;
                }
              }
            }
          }
        }
        
        #Links that are not to Uniprot will use double time he Head function and finally the Net::FTP check
        #NOT UNIPROT
         if (!head($configHash->{'enzymeDBLink'}) ) {
          if (!head($configHash->{'enzymeDBLink'})) {
            if (testFTP($configHash->{'enzymeDBLink'})!= 1) {
              annoPrint ("Please check if [".$configHash->{'enzymeDBLink'}."] is a correct URL. Annocript will continue...\n");#
              $wrongLinks++;
            }	
          }
        }
        if (!head($configHash->{'cdDBLink'})) {
          if (!head($configHash->{'cdDBLink'})) {
            if (testFTP($configHash->{'cdDBLink'}) != 1) {
              annoPrint ("Please check if [".$configHash->{'cdDBLink'}."] is a correct URL. Annocript will continue...\n"); #
              $wrongLinks++;
            }
          }
        }
        #This link always has some more troubles... let's add one more check!
        if (!head($configHash->{'cdTableLink'})) {
          if (!head($configHash->{'cdTableLink'})) {
            if (!head($configHash->{'cdTableLink'})) {
              if (testFTP($configHash->{'cdTableLink'}) != 1) {
                annoPrint ("Please check if [".$configHash->{'cdTableLink'}."] is a correct URL. Annocript will continue...\n"); #
                $wrongLinks++;
              }
            }
          }
        }
        
        if (!head($configHash->{'rfamDBLink'})) {
          if (!head($configHash->{'rfamDBLink'})) {
            if (testFTP($configHash->{'rfamDBLink'}) != 1) {
              annoPrint ("Please check if [".$configHash->{'rfamDBLink'}."] is a correct URL. Annocript will continue...\n"); 
              $wrongLinks++;
            }
          }
        }

        if (!head($configHash->{'GOTermsLink'})) {
          if (!head($configHash->{'GOTermsLink'})) {
            if (testFTP($configHash->{'GOTermsLink'}) != 1 ) {
                annoPrint ("Please check if [".$configHash->{'GOTermsLink'}."] is a correct URL. Annocript will continue...\n"); #
                $wrongLinks++;
            }
          }
        }

        if (!head($configHash->{'pathwaysTableLink'})) {
          if (!head($configHash->{'pathwaysTableLink'})) {
            if ( testFTP($configHash->{'pathwaysTableLink'}) != 1 ) {
              annoPrint ("Please check if [".$configHash->{'pathwaysTableLink'}."] is a correct URL. Annocript will continue...\n"); #
              $wrongLinks++;
            }
          }
        }
       
      if ($wrongLinks > 0){
        my $text = "Since there is some problem with connection, Annocript is trying different kind of download. Be sure that the links above are working".
        " and then press y to continue.\n Continue? (y or n)\n";
        #my $text ="The links above seems to not work. You should check them in a browser and if (1) you are sure they work or (2) you never".
        #" modified this links in the configuration file, then there is something in your connection that is not working. In this last case".
        #" since Annocript uses more than one technique to download the files you can leave it running. Continue? (y or n)\n"; 
        my $choose = get_input($text,'^[yYnN]$');
        annoPrint($choose."\n",1);#Print on the log file
          if ( ($choose eq "n") or ($choose eq "N")){ 
           die "Ok. Let's correct the links and restart Annocript! Exiting...\n";
          }        
       }  
      annoPrint (". Everything is ok!\n");
      
      		     
      
  }#This is a simple control that is done even if the database is not created. It is more simple and faster.
  elsif (! (checkLink($configHash->{'rfamDBLink'}) and checkLink($configHash->{'cdDBLink'}) and checkLink($configHash->{'enzymeDBLink'}) and 
          checkLink($configHash->{'GOTermsLink'}) and checkLink($configHash->{'GODBLink'})
            and checkLink($configHash->{'unirefVerLink'}) and checkLink($configHash->{'uniprotVerLink'}) and checkLink($configHash->{'uniprotVerLink'}) 
              and checkLink($configHash->{'tremblDBLink'}) and checkLink($configHash->{'swissprotDBLink'}) ) ){
                annoPrint ("Pay attention... you did some error in writing links in the configuration file. Please check it...\n");
  }elsif ( ($configHash->{'doBuildOutput'} eq 'YES') or ($configHash->{'doBuildOutput'} eq 'YES' )){
    
	  #The first check is on the permissions. If the user does not have, he cannot create database or perform analysis
	  if (!check_DB_user_permissions($configHash->{'platform'},$configHash->{'host'},
				  $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'})) {
		 annoDie("ERROR: You do not have enough privileges to create the DB. Please ask your adiministrator. Annocript will close...\n");
	  }
  }
  
  
  #Check if the database chose for the mostly expressed domains is correct
  my $cdName4Expression = $configHash->{'cdName4Expression'};
  #print "cdName4Expression : $cdName4Expression\n";
  if ( !(grep {/$cdName4Expression/ } @cdDBAdmitted) ){
    annoDie("Can't find CD database $cdName4Expression. You should "
        ."choose one among:\n ".print_array(\@cdDBAdmitted)."\nAnnocript will exit...\n");
  }
      
  #Check on the log file for 
  if ( !($configHash->{'headParsNAValues'} =~ /\.log?/)){
    annoDie( $configHash->{'headParsNAValues'}." is not a correct log file name. It must be something like: file.log\n");
  }
 
  if ( !($configHash->{'uniprotGenesNotPresent'} =~ /\.log?/)){
    annoDie( $configHash->{'uniprotGenesNotPresent'}." is not a correct log file name. It must be something like: file.log\n");
  }
  
  #A check on the user-id and password: it can be of every type of character so 
  #it don't needs a particular control
  if ( !(defined ($configHash->{'mySqlUser'})) or !(defined ($configHash->{'mySqlPass'}) ) ){
    annoDie( "ERROR: You need to set an user-id and password for MySQL. Please contact your adiministrator".
    " if you can't. Annocript will exit...\n");
  }#else{annoPrint ("MYSQL account: ".$configHash->{'mySqlUser'}." and pwd: ".$configHash->{'mySqlPass'}."\n";}#DEBUGCODE
  
  #A check on the parsing file path
  if ( -e $dataFolder.$configHash->{'parsingFilePath'} ){
    annoDie( "ERROR : File ".$configHash->{'parsingFilePath'}." already exists.\n"); #Please change the parsingFilePath in ".$configFilePath."\n";
  }
              

  #If the user sets the config file to not execute the DB_CREATION, this means that nothing will 
  #be downloaded: the user has already a database.
  my $dbDataFolder; 
    if ( $configHash->{'doDbCreation'} ne "YES"){
        if( !(scalar <$mainDataFolder/*>) ){
           annoDie( "ERROR: data folder is empty... You must run Annocript with doDbCreation = YES. Closing");
        }else{
          #In the first case Annocript will ask to the user what folder to use as database
          annoPrint ("DB Creation not active: Database will not be downloaded and installed.\n");
          my $useDb = get_input("You chose to use an existent database, are you sure?(y or n): ",'^[YynN]$');
          annoPrint($useDb."\n",1);
          if ( ($useDb eq "y") or ($useDb eq 'Y') ) {
						##User wants to choose another folder
						#Opens the dir with all the files
						opendir DIR, $mainDataFolder or annoDie( "ERROR : cannot open dir $mainDataFolder");
						my @files = readdir DIR;#All the files are in this list  
						print "Available database versions:\n";
									foreach my $file (@files){
							if ($file =~ /\d_\d/){
								print $file." ";
							}
						}
						print "\n";
									annoPrint ("Choose the database version: ");
									my $myFolder = <STDIN>;
									chomp $myFolder;
									while( !(-d $mainDataFolder."/".$myFolder) ){
										foreach my $file (@files){
							if ($file =~ /\d_\d/){
								print $file." ";
							}
							
							}
							print "\n";
										annoPrint ("Choose the database version: ");
										$myFolder = <STDIN>;
										chomp $myFolder;
									}
									$dbDataFolder = $mainDataFolder."/".$myFolder;
									annoPrint($myFolder."\n",1);   
									
								 
									#Once that the folder has been chose we need to ask to the user if he wants
									#to use UNIREF or UNIPROT_KB databases. 
									my $uniprotDBName = "uniprotkb_".$myFolder;
									my $unirefDBName = "uniref_".$myFolder; 
									
									#Here we access the mysql system to see if the databases are present both or uniquely
									my $uniprotDBPresent = db_present($uniprotDBName,$configHash->{'platform'},$configHash->{'host'},
									$configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
									my $unirefDBPresent =  db_present($unirefDBName,$configHash->{'platform'},$configHash->{'host'},
									$configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
									
									#If they both don't exist
									if( ( $uniprotDBPresent == 0)  and ($unirefDBPresent == 0) ){
									 annoDie( "No database in folder $dbDataFolder... You must run Annocript with doDbCreation = YES. Closing");
									#They both are present
									}elsif (( $uniprotDBPresent > 0)  and ($unirefDBPresent > 0) ){
										my $text ="UniRef and TrEMBL databases are both there. What do you want to use? UniRef (u) or TrEMBL (t)?\n".
										"(Suggestion: if you are using a file with organisms to match with, you can only use TrEMBL database.)\n"; 
										my $choose = get_input($text,'^[uUtT]$');
										annoPrint($choose."\n",1);#Print on the log file
										if ( ($choose eq "u") or ($choose eq "U")){ 
											$configHash->{'dbInUse'}= "uniref";
											$configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
											 $configHash->{'unirefDB'} = choose_uniref_database($dbDataFolder);#extract_name($configHash->{'unirefDBLink'},'gz');
										}else{ 
											$configHash->{'dbInUse'}= "uniprotkb";
											$configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
											
										}
										#$configHash->{'swissProtDB'} = extract_name($configHash->{'swissprotDBLink'},'gz');
									}#Only uniprot kb is present
									 elsif ( ($uniprotDBPresent > 0)  and ($unirefDBPresent == 0) ){
											annoPrint ("UniProt database is the only present. Using it...\n");
											$configHash->{'dbInUse'}= "uniprotkb";
											$configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
							
										#Only uniref is there
										}elsif ( ($uniprotDBPresent == 0)  and ($unirefDBPresent > 0) ){
											annoPrint ("UniRef database is the only present. Using it...\n");
											$configHash->{'dbInUse'}= "uniref";
											 $configHash->{'unirefDB'} = choose_uniref_database($dbDataFolder);#extract_name($configHash->{'unirefDBLink'},'gz');
											$configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
										}
										
									#We store the names for all the database. They will not be downloaded but we need them later  
								 
									$configHash->{'tremblDB'} = extract_name($configHash->{'tremblDBLink'},'gz');
									$configHash->{'swissProtDB'} = extract_name($configHash->{'swissprotDBLink'},'gz');
									
									$configHash->{'database'} = $configHash->{'dbInUse'}.'_'.$myFolder;#Gives the name to the database
									$configHash->{'uniprotKBDB'} = "uniprot_kb_".$myFolder.".fasta";#Creates the name for the uniprot kb
										 
									$configHash->{'cdDB'} = choose_CD_database($dbDataFolder);#The user choose the CD database if it is present
									$configHash->{'cdDBPath'} = $dbDataFolder."/".$configHash->{'cdDB'};
									$configHash->{'dbVersion'} = $myFolder;	
          }
          #The user has to change the config file	
          else { annoDie( "Please change parameter doDbCreation in the config file. Annocript will close...\n");}	
     
     #Then do many checks on the indexed databases and on the choose of a single or many organisms but only
     #if the user wants to execute analyses, otherwise it does not matter
     if ($configHash->{'doExecutePrograms'} eq 'YES'){     
			 #Changing the database name if a single organims is to be used, only when using Trembl
			 if ( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniprotkb") ){
				 #If the file with the sequences exists in the working folder... 
				if (-e $workingFolder.'/'.$configHash->{'blastedOrganism'}){
					my $blOrg = extract_name($configHash->{'blastedOrganism'},0);
					
					$configHash->{'tremblDB'} = "TR_$blOrg.fasta";
					$configHash->{'swissProtDB'} = "SP_$blOrg.fasta";
					#check if there is a file named in same way in the folder. if it is there, the user is alerted
					if (-e $dbDataFolder."/".$configHash->{'tremblDB'}){
						my $sure= get_input("The indexed database named ".$configHash->{'tremblDB'}." is present in".
						" $dbDataFolder. Shall I use it?(y or n)",'^[YynN]$');
						annoPrint($sure."\n",1);
							if ( ( $sure eq "n") or ( $sure eq "N")){
							annoPrint(	"Ok, please change the name in the configuration file...\n");
							exit;
							}else{annoPrint(	"Ok, Annocript will use it or overwrite if not complete...\n");}
					}
					if (-e $dbDataFolder."/".$configHash->{'swissProtDB'}){
						my $sure= get_input("The indexed database named ".$configHash->{'swissProtDB'}." is present in".
						" $dbDataFolder. Shall I use it?(y or n)",'^[YynN]$');
						annoPrint($sure."\n",1);
							if ( ( $sure eq "n") or ( $sure eq "N")){
							annoPrint(	"Ok, please change the name in the configuration file...\n");
							exit;
							}else{annoPrint(	"Ok, Annocript will use it or overwrite if not complete...\n");}
							#If for some reason the database is not complete, you will overwrite it
					}
					copy($workingFolder.'/'.$configHash->{'blastedOrganism'},$dbDataFolder.'/'.$configHash->{'blastedOrganism'});				
				}else {
					annoDie("There is no file named ".$configHash->{'blastedOrganism'}." in $workingFolder. Please".
					" check you wrote a correct name in the config file and that the file with organisms names is in your".
					" working folder. If you want to blast against all organisms please change the value with 'all'. \n");
					}  
				}elsif( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniref")){
					annoPrint("When using UniRef all the organisms will be evalutated for the blast. If you want to blast against ".
					" specific organisms you should use TrEMBL.\n");
						my $text ="Continue? (y or n)"; 
					my $choose = get_input($text,'^[yYnN]$');
					annoPrint($choose."\n",1);#Print on the log file
					#The user choose to continue using UniRef
					if ( ($choose eq 'n') or ($choose eq 'N')) {
						annoPrint("Please change parameter 'blastedOrganism' or choose Trembl when asked. Annocript will stop...\n");
						exit;
					}
				 }
				
				#Check on the databases. The indexed should be present, otherwise the user will be alerted that
				#he should start dbCreation procedure
				if ($configHash->{'dbInUse'} eq "uniprotkb"){
					if ( (indexed_db_present($configHash->{'tremblDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
					 annoDie("Database ".$configHash->{'tremblDB'}." is not indexed. Please restart annocript with doDbCreation=YES\n");
					}				
				}
				
				if ($configHash->{'dbInUse'} eq "uniref"){
					if ( (indexed_db_present($configHash->{'unirefDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
					 #annoPrint("Database ".$configHash->{'uniref'}." is not indexed. Annocript will do it!\n");
					 annoDie("Database ".$configHash->{'unirefDB'}." is not indexed. Please restart annocript with doDbCreation=YES\n");
					}				
				}
				
				if ((indexed_db_present($configHash->{'swissProtDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
				 annoDie("Database ".$configHash->{'swissProtDB'}." is not indexed. Please restart annocript with doDbCreation=YES\n");
				}
		 }
     }
  }else{
    #The database has to be created. The user can choose what (UNIREF or UNIPROT_KB)
    my $text ="You chose to create the database. Do you want to create with UniRef (u) or TrEMBL (t)? (u or t)\n".
    "[Suggestion: if you want to blast against specific organisms you should use TrEMBL]\n"; 
    my $choose = get_input($text,'^[uUtT]$');
    annoPrint($choose."\n",1);#Print on the log file
    #The user choose to create uniref
    if ( ($choose eq "u") or ($choose eq "U")){ 
      $configHash->{'dbInUse'}= "uniref";
      $configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
      $configHash->{'unirefDB'} = extract_name($configHash->{'unirefDBLink'},'gz');
    }
    #The user choose to create uniprot_kb
    else{ 
      $configHash->{'dbInUse'}= "uniprotkb";
      $configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
      $configHash->{'tremblDB'} = extract_name($configHash->{'tremblDBLink'},'gz');
    }
    #The name of the SwissProt db
    $configHash->{'swissProtDB'} = extract_name($configHash->{'swissprotDBLink'},'gz');
    
 
		  
    #To reinstall an already created db uncomment the first and comment the second following lines
    #my $dbVersion = '2013_08';
    
    my $dbVersion = getDBVersion($configHash->{'versionDBLink'},$mainDataFolder,$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'}, $configHash->{'dbInUse'});
    annoPrint ("\ndbVersion: ".$dbVersion."\n");
  	$configHash->{'dbVersion'} = $dbVersion;
    		
    #..and the folder.
    $dbDataFolder = $mainDataFolder."/".$dbVersion;#Create the name for the folder..
    $configHash->{'database'} = $configHash->{'dbInUse'}.'_'.$dbVersion;#Give the name to the database
    annoPrint ("Folder for DB will be: ".$dbDataFolder." and database name: ".$configHash->{'database'}."\n");
    $tempTableFolder =  $dbDataFolder."/".$tempTableFolder;#Creating the tables folder
    $configHash->{'tempTableFolder'}  = $tempTableFolder;#Creating a global variable to share among modules
    
    #Creating the db folder
    unless(-d $dbDataFolder){
      annoPrint ($dbDataFolder." doesn't exists. Creating folder $dbDataFolder...\n");
      mkdir $dbDataFolder or annoDie( "ERROR : Can't create $dbDataFolder. Check permissions!\n");
    }else{#If the folder is there, check the UniRef db he wants to install
	  if ($configHash->{'dbInUse'} eq "uniref"){
		  #Check if there is a uniref database in the db folder, if something is there the version of uniref is compared
		  #with the uniref link and if it is different the user is alerted
		  my @ufPresent = check_uniref_presence($dbDataFolder);
		  if (scalar(@ufPresent) > 0){
			  my $fromLink = extract_name($configHash->{'unirefDB'},1);
			  if( !(grep  {/$fromLink/} $ufPresent[0]) ){
				annoDie( "The type of UniRef you want to download is not the same as the one you downloaded in a previous execution".
				" of Annocript. Annocript can manage only one type of UniRef (50,90 or 100) for each versions of the Uniprot db.".
				"Check the unirefDBLink and unirefVerLink in CONFIGURATION/config_annocript.txt file.\n");
			  }
		  }
	  }
	  annoPrint ("A folder $dbDataFolder is already existent. Annocript will continue...\n ");
	}
	
 #Creates a folder where all the tables will be put
	unless(-d $tempTableFolder){
		annoPrint ($tempTableFolder." doesn't exists. Creating folder $tempTableFolder...\n");
		mkdir $tempTableFolder or annoDie( "ERROR : Can't create $tempTableFolder. Check permissions!\n");
	}else{annoPrint ("A folder $tempTableFolder is already existent. Annocript will continue...\n ");}
	
	$configHash->{'uniprotKBDB'} = "uniprot_kb_".$dbVersion.".fasta";#Create the name for the uniprot kb
	
	#Move the file with the database version to the folder just created
	move($mainDataFolder."/".(extract_name($configHash->{'versionDBLink'},0)),$dbDataFolder); 
		
	#Give the  name to the CDD database extracting from the link
	$configHash->{'cdDB'} = extract_name($configHash->{'cdDBLink'},"1");
	#Sometimes in uniprot they put the name of the zip file with LE (Little Endian) and BE (Big Endian)
	$configHash->{'cdDB'} = clean_cd_name($configHash->{'cdDB'});
	 
	$configHash->{'cdDBPath'} = "$dbDataFolder/".$configHash->{'cdDB'};
	
	annoPrint ("Conserved domains database: ".$configHash->{'cdDB'}."\n");
	$text ="Annocript will download several files for the database generation and BLAST execution. ".
	" Do you want to remove the file that Annocript will not use for the analysis? (y or n) [Suggested: n].\n"; 
	$choose = get_input($text,'^[YynN]$');
	annoPrint( $choose."\n",1);
	if ( ($choose eq "y") or ($choose eq "Y")){ 
		$configHash->{'eraseFiles'}= "YES";
	}
	elsif ( ($choose eq "n") or ($choose eq "N")) { 
		$configHash->{'eraseFiles'}= "NO";
	}else {annoDie( "ERROR : Exiting Annocript...\n");}
	annoPrint($choose."\n",1);
	
	#Here Annocript decides if to start the DB creation using more memory and doing fast or not
	#It depends on how much memory the system has but it can be choose inside the annocript configuration file
	my $ram = ram_memory();
	if ( $ram/1024/1024 > $configHash->{'bigMemory'} ){
		$configHash->{'useHugeRAM'} = 'YES';
		annoPrint("You have enough memory to create more rapidly the db. Annocript will use a lot of memory now! :-D\n");
	}else{
		$configHash->{'useHugeRAM'} = 'NO';
		annoPrint("You don't have enough memory to be more rapid. Annocript will be little slow but won't use lot of memory! :-)\n");
	}
      
     #Changing the database name if a single organims is to be used
	 if ( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniprotkb") ){
	  #Search the file with organisms	
	  if (-e $workingFolder.'/'.$configHash->{'blastedOrganism'}){
		my $blOrg = extract_name($configHash->{'blastedOrganism'},0);
		$configHash->{'tremblDB'} = "TR_$blOrg.fasta";
		$configHash->{'swissProtDB'} = "SP_$blOrg.fasta";
		#check if there is a file named in same way in the folder. if it is there, the user is alerted		
	    if (-e $dbDataFolder."/".$configHash->{'tremblDB'}){
			annoPrint("An existing database ".$configHash->{'tremblDB'}." is present in the folder ".
			"$dbDataFolder. Please change the name of the file of organisms or remove this file.\n"  );
			exit;
		}
		if (-e $dbDataFolder."/".$configHash->{'swissProtDB'}){
			annoPrint("An existing database ".$configHash->{'swissProtDB'}." is present in the folder ".
			"$dbDataFolder. Please change the name of the file of organisms or remove this file.\n"  );
			exit;
		}		
		#Shift the file with the organisms to blast inside the db folder
		copy($workingFolder.'/'.$configHash->{'blastedOrganism'},$dbDataFolder.'/'.$configHash->{'blastedOrganism'});
				
	  }else {
		annoDie("There is no file named ".$configHash->{'blastedOrganism'}." in $workingFolder. Please".
		 " check you wrote a correct name in the config file and that the file with organisms names is in your".
		 " working folder. \n");
	  } 
	 }elsif( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniref")){
			annoPrint("When using UniRef all the organisms will be evaluated for the blast. If you want to blast against ".
				" specific organisms you should use Trembl.\n");
			my $text ="Continue? (y or n)"; 
			my $choose = get_input($text,'^[yYnN]$');
			annoPrint($choose."\n",1);#Print on the log file
			#The user choose to continue using UniRef
			if ( ($choose eq 'n') or ($choose eq 'N')) {
				annoPrint("Please change parameter 'blastedOrganism' or choose Trembl when asked. Annocript will stop...\n");
				exit;
			}
	  }
			    
	#}
  }#END OF CODE WHEN THE DATABASE IS TO BE GENERATED
  	
  if (defined $dbDataFolder){
    $configHash->{'dbDataFolder'} = $dbDataFolder;
  }
  
  #A check on the database name...
  uniprot_db_selection( $configHash->{'database'},$configHash->{'platform'},$configHash->{'host'},
  $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'},
  $configHash->{'doDbCreation'},$configHash->{'doBuildOutput'});
}



=head2 check_DB_user_permissions

 Title   : check_DB_user_permissions
 Usage   : check_DB_user_permissions(  - local folder = the folder where actually is the GFF3_AND_OUTPUT script
                               );

 Function:  this function accesses the database with the user and password given and says if the user has sufficient privileges to 
        create a database
 
 Returns : 1 if check succeed -1 otherwise

=cut
sub check_DB_user_permissions{
	my $platform = shift;
	my $host = shift;
	my $port = shift;
	my $mySqlUser = shift;
	my $mySqlPass = shift;
  
  
  my $retVal = -1;
  my $privilegesType = 'GRANT ALL PRIVILEGES ON';
  
  annoPrint("Checking your privileges for MySQL...");
  # PERL DBI CONNECT
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } ) 
    or die "Unable to connect (Check mySqlUser and mySqlPass in the config file): $DBI::errstr\n";;
  
  #Query to show grants for the user
  my $query = "SHOW GRANTS FOR '$mySqlUser'\@'$host'";
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my  @data = $sth->fetchrow_array();
  
  my $privileges = $data[0];
  #annoPrint( "$privileges;\n\n");
   
  if($privileges =~ /$privilegesType/){
    $retVal = 1;
    print "DONE!\n";
  }
}


=head2 choose_CD_database

 Title   : choose_CD_database
 Usage   : choose_CD_database(  - local folder = the folder where actually is the GFF3_AND_OUTPUT script
                               );

 Function:  this function searches in the data folder for the CD databases present. If something is there then it will show
          and the user can select one of them. To check if a CD database exists it controls if a set of extensions is present
          If not the database will not be shown in the list. 
 
 Returns : The name of the database chose or Null if nothing is there

=cut
sub choose_CD_database{
  my $dbDataFolder = shift;
  
  
  my $cdDBchose = '';#Return value. A db or null
  my $cdDBFound = 0;#If a db is found is 1
  my @cdDBs = ();#cdDBs found
  
  #Check if the CD db exists
  #Opens the dir with all the databases
  opendir DIR, $dbDataFolder or annoDie( "ERROR : cannot open dir $dbDataFolder");
  my @files= readdir DIR;#All the files are in this list
  
  #A control starts on each file of the folder        
  foreach my $file (@files){
    #If a .rps file is found then it will search for the other extensions
    my $ext1 = $cdFilesExts[0];#This is .rps
    if ($file =~ /$ext1/){
      $file = (split(/\./,$file))[0];
      # print "$file matches $ext1\n";
      my $matched = 0;
      #If for each needed extension the file is there
      foreach my $ext (@cdFilesExts){
        my $fToCheck = $dbDataFolder.'/'.$file.$ext;
        #print "Check if $fToCheck exists\n";
        if (-e $fToCheck){
          #print "Wow.. exists: ".$dbDataFolder.'/'.$file.$ext."\n";
           $matched++;
        }
      }
     #IF the matches are as much as the extensions  
     if ($matched == scalar(@cdFilesExts)){
       #print "WOW... all files are there!!\n";#DEBUGCODE
       #Keep a list with all the CDDBs
       push(@cdDBs, $file);
       #and print on screen the file to choose
       $cdDBFound = 1;
      }
    }
  }
 #If nothing is there Annocript should die       
  if ($cdDBFound == 0) {
    annoDie ("You choosed to use existent conserved domains db but nothing has been found in $dbDataFolder.".
        " Please restart Annocript with dbCreation = YES. Annocript will go out for other commitments...\n"); 
  }else{
    #only one database is present
    if (scalar(@cdDBs) == 1){
      print $cdDBs[0]." is the only CD database present. Using it...\n";
      $cdDBchose = $cdDBs[0];
    }else{
      #Print the cd database names until the user don't choose a correct one
       annoPrint ("Choose your cd database: ");
       print_array(\@cdDBs);
       $cdDBchose = <STDIN>;
       chomp $cdDBchose;
       while( !( grep {/\b$cdDBchose\b/} @cdDBs) ){
          annoPrint ("Choose your cd database: ");
          print_array(\@cdDBs);
          $cdDBchose = <STDIN>;
          chomp $cdDBchose;
       }
       annoPrint ("Cd database chose:  $cdDBchose. Here we go!\n");
      }
  }
  
  return $cdDBchose;
}

=head2 check_uniref_presence

 Title   : check_uniref_presence
 Usage   : check_uniref_presence( 
							- dbDataFolder = the folder where actually is to find the db
                               )

 Function: this function searches in the data folder for the UniRef databases present. If something is there then it will 
 return
 
 Returns : The name of the database present

=cut
sub check_uniref_presence{
  my $dbDataFolder = shift;
  
  my $uniref50 ='uniref50.fasta';
  my $uniref90 ='uniref90.fasta';
  my $uniref100 ='uniref100.fasta';
  my $u50ok = 0;
  my $u90ok = 0;
  my $u100ok = 0;
  my @uDBpresent = ();
    
  #Create an array with versions
  my @ufVersions = ('50','90','100');
  
  #Check if the  db exists
  #Opens the dir with all the databases
  opendir DIR, $dbDataFolder or annoDie( "ERROR : cannot open dir $dbDataFolder");
  my @files = readdir DIR;#All the files are in this list
  
  #A control starts on each file of the folder        
  foreach my $file (@files){
    #If a .rps file is found then it will search for the other extensions
    #my $ext1 = $ufVersions[0];#This is .rps
    
    if ($file eq $uniref50){
		$u50ok = 1;
		push(@uDBpresent,$uniref50);
	}
    if ($file eq $uniref90){
		$u90ok = 1;
		push(@uDBpresent,$uniref90);
	}
    if ($file eq $uniref100){
		$u100ok = 1;
		push(@uDBpresent,$uniref100);
	}
  }
  
  return @uDBpresent;
}



=head2 choose_uniref_database

 Title   : choose_uniref_database
 Usage   : choose_uniref_database( - local folder = the folder where actually is the GFF3_AND_OUTPUT script
                               );

 Function: this function searches in the data folder for the UniRef databases present. If something is there then it will show
          and the user can select one of them. To check if a CD database exists it controls if a set of extensions is present
          If not the database will not be shown in the list. 
          Actually this function searches the UniRef database type you have but Annocript does not permit more than one type
          of UniRef in the same folder. I leave this function as it is so that if in the future we need to design
          a new advanced Annocript with possibility to use more UniRef types it will work fine also on that!
 
 Returns : The name of the database chosen or Null if nothing is there

=cut
sub choose_uniref_database{
  my $dbDataFolder = shift;
   
  my $dbchosen = '';#Return value. A db or null
  my $dbFound = 0;#If a db is found is 1
  my @dbs = ();#dbs found
  
  my $uniref50 ='uniref50.fasta';
  my $uniref90 ='uniref90.fasta';
  my $uniref100 ='uniref100.fasta';
  my $u50ok = 0;
  my $u90ok = 0;
  my $u100ok = 0;
  my @uDBpresent = ();
    
  #Create an array with versions
  my @ufVersions = ('50','90','100');
  
  #Check if the  db exists
  #Opens the dir with all the databases
  opendir DIR, $dbDataFolder or annoDie( "ERROR : cannot open dir $dbDataFolder");
  my @files = readdir DIR;#All the files are in this list
  
  #A control starts on each file of the folder        
  foreach my $file (@files){
    #If a .rps file is found then it will search for the other extensions
    #my $ext1 = $ufVersions[0];#This is .rps

    
    if ($file eq $uniref50){
		$u50ok = 1;
		push(@uDBpresent,$uniref50);
	}
    if ($file eq $uniref90){
		$u90ok = 1;
		push(@uDBpresent,$uniref90);
	}
    if ($file eq $uniref100){
		$u100ok = 1;
		push(@uDBpresent,$uniref100);
	}
  }
  
 #If nothing is there Annocript should die       
  if ( scalar(@uDBpresent) == 0 ) {
    annoDie ("No uniref database present. Please restart Annocript with dbCreation = YES. Exiting...\n"); 
  }else{
    #only one database is present
    if ( scalar(@uDBpresent) == 1){
      print $uDBpresent[0]." is the only uniref database present. Using it...\n";
      $dbchosen = $uDBpresent[0];
    }else{
      #Print the cd database names until the user doesn't choose a correct one
       annoPrint ("Choose your UniRef database: ");
       print_array(\@uDBpresent);
       $dbchosen = <STDIN>;
       chomp $dbchosen;
       while( !( grep {/\b$dbchosen\b/} @uDBpresent) ){
          annoPrint ("Choose your cd database: ");
          print_array(\@uDBpresent);
          $dbchosen = <STDIN>;
          chomp $dbchosen;
       }
       annoPrint ("UniRef database chosen:  $dbchosen. Here we go!\n");
      }
  }
  
  return $dbchosen;
}

=head2 checkGFF3_AND_OUTPUT

 Title   : checkGFF3_AND_OUTPUT
 Usage   : checkGFF3_AND_OUTPUT(  
								- local folder = the folder where actually is the GFF3_AND_OUTPUT script
                               );

 Function: Checks the config hash parameters regarding the ouput creation: how they are written, if they are present and so on..
 Returns : nothing

=cut
sub checkGFF3_AND_OUTPUT {
  my $localFolder = shift;
  
  #A check on the GFF database ...
  my $name = $configHash->{'name'};
  #Length max to get from a file to give the name to the GFF db
  if ( length($name) > 15 ){
      $name = substr($name,0,15);
   }
  
  #Remove strange chars from the name and subs
  $name =~ s/[\-\(\)\s\_\/\.]/_/g;
  $configHash->{'gffDB'} = $configHash->{'gffDB'}.'_'.$name;
  
  #This is a control on the presence of the GFF database. If it exists means that the user previously selected to not destroy it.
  detect_gff_db( $configHash->{'gffDB'},$configHash->{'platform'},$configHash->{'host'},
              $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
		
  #A check of the file extensions...
  if ( !($configHash->{'blastxGFF3FileName'} =~ /\.gff?/)){
    annoDie($configHash->{'blastxGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n");
  }  
  if ( !($configHash->{'blastxTRGFF3FileName'} =~ /\.gff?/)){
    annoDie($configHash->{'blastxTRGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n");
  } 
  if ( !($configHash->{'blastxSPGFF3FileName'} =~ /\.gff?/)){
    annoDie( $configHash->{'blastxSPGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n");
  } 
  if ( !($configHash->{'rpstblastnGFF3FileName'} =~ /\.gff?/)){
    annoDie( $configHash->{'rpstblastnGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n");
  } 
  if ( !($configHash->{'blastnGFF3FileName'} =~ /\.gff?/)){
    annoDie( $configHash->{'blastnGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n");
  }
    
  if ( !($configHash->{'ORFFastaFileName'} =~ /(\.fasta?)|(\.fa?)/)){
    annoDie( $configHash->{'ORFFastaFileName'}." is not a correct FASTA file name. It must be something like: file.fasta\n");
  } 
  if ( !($configHash->{'fastaSeqs'} =~ /(\.fasta?)|(\.fa?)/)){
    annoDie( $configHash->{'fastaSeqs'}." is not a correct FASTA file name. It must be something like: file.fasta\n");
  }
  if ( !($configHash->{'outFileName'} =~ /\.txt?/)){
    annoDie( $configHash->{'outFileName'}." is not a correct txt file name. We prefer to have a txt file for all results: output.txt\n");
  }
  #A control on the type of the output of blast programs
  if ( !($configHash->{'blastOutType'} =~ /blast/)){
    annoDie( $configHash->{'blastOutType'}." the blast out type supported is blast, keep in mind and modify the configuration file!\n");
  }
  
  #Checking the variable goTermsAss. It can have only one of three values
  if ( !($configHash->{'goTermsAss'} =~ /(proteins?)|(domains?)|(both?)/)){
    annoDie( "You can choose only one between [proteins, domains] as a value for goTermAss in the configuration file. ".$configHash->{'goTermsAss'}." is not allowed!\n");
  } 
}

=head2 check_programs_presence

 Title   : check_programs_presence
 Usage   : check_programs_presence(  - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function: Execute a set of controls on the presence of programs depending from the user choices
          done in the configuration file.
 Returns : nothing

=cut
sub check_programs_presence {
	
  #Here we collect the path to the programs using an array initialized in the beginning of the page
  # with all the programs used by Annocript
  foreach my $program (@blastProgramsInUse){
    my $var= $program.'Path';
    $configHash->{$var} = $configHash->{'blastPath'}.$program;
  }
  
	    #CONTROLS THE PATH TO THE PROGRAMS
    if (!(check_presence($configHash->{'rpstblastnPath'}))){
      annoDie($configHash->{'rpstblastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    } 
    
    if (!(check_presence($configHash->{'rpsblastPath'}))){
      annoDie($configHash->{'rpsblastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
    
    if ( (!check_presence($configHash->{'blastxPath'})) ){
      annoDie($configHash->{'blastxPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");  
    } 
    
    if ( (!check_presence($configHash->{'blastpPath'})) ){
      annoDie($configHash->{'blastxPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");  
    } 
    
    if ( (!check_presence($configHash->{'blastnPath'})) ){
      annoDie($configHash->{'blastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
    
    if ( (!check_presence($configHash->{'tblastnPath'})) ){
      annoDie($configHash->{'blastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
    
    if ( (!check_presence($configHash->{'makeblastdbPath'})) ){
      annoDie($configHash->{'makeblastdbPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
    
    if ( (!check_presence($configHash->{'portraitPath'})) ){
      annoDie($configHash->{'portraitPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
    
    if ( (!check_presence($configHash->{'dna2pepPath'}))  ) {
      annoDie($configHash->{'dna2pepPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n");
    }
}

=head2 check_dependencies

 Title   : check_dependencies
 Usage   : check_dependencies(  - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function: Execute a set of controls on the output and presence of programs depending from the user choices
          done in the configuration file.
 Returns : nothing

=cut
sub check_dependencies {
  #If the programs have to be executed paths are controlled and is checked if some outputs are present
  if ($configHash->{'doExecutePrograms'} eq 'YES'){  
    #We check if in the blast out folders there is some BLAST output
    if( (scalar <$outFolder/rpstblastn/*>) or (scalar <$outFolder/blastx/*>) or (scalar <$outFolder/blastn/*>) ) {
      my $removeOut = get_input("\n WARNING: Some output files already exist! If you will run a program (BLAST, PORTRAIT, DNA2PEP) the corresponding ".
                                "output will be overwritten. Instead, output of BLASTx against Uniref will be different from out of BLASTx against Uniprot_kb.".
                                " Continue?(y or n)",'^[YynN]$');
      annoPrint($removeOut."\n",1);#Print on the log file
      if (($removeOut eq "n") or ($removeOut eq "N") ){
            annoDie( "Go in the output folder and save or cancel the output you want/don't. Or choose programs to be executed. Exiting...\n ");
      }else { print	"Ok, Annocript will continue...\n";}
    }
  #If programs will not be run and the user wants to print the final output table then a check is needed if outputs of programs
  #are present
  }elsif ( ($configHash->{'doExecutePrograms'} ne 'YES') and ($configHash->{'doBuildOutput'} eq 'YES') ){
    if( !(scalar <$outFolder/rpstblastn/*>) and !(scalar <$outFolder/blastx/*>) and !(scalar <$outFolder/blastn/*>)) {
      annoDie( "Blast output folders are empty... Please run the programs to obtain annotations. Exiting...\n");
    }
  #If only db creation is active we only check the presence of the makeblastdb program
  }elsif ( $configHash->{'doDbCreation'} eq 'YES'){
    if (!check_presence($configHash->{'makeblastdbPath'})){
        annoDie("You miss makeblastdb. Please check the links in the configuration file or re-install the BLAST program. Exiting...\n");
    }
  } 
  
  #If nothing is there in the folders we cannot start the output and stats creation
  if ( ($configHash->{'extractStatistics'} eq 'YES') and (($configHash->{'doBuildOutput'} ne 'YES')) ){
            if (!check_presence($outFolder.'/'.$configHash->{'outFiltered'})){
                annoDie("You cannot start the statistics without the filtered ouput of Annocript. Please restart with doBuildOutput = YES\n");
            }
  }
  
}


=head2 check_Annocript_integrity

 Title   : check_Annocript_integrity
 Usage   : check_Annocript_integrity( - configHash = the piece of the config hash for PROGRAMS_EXEC
							 - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function:  checks variables mostly inside the annocript configuration file. 
 Returns : nothing

=cut
sub check_Annocript_integrity {
    
  #######CHECKING integrity of Annocript
  if ( ! ( (-e $configHash->{'usefulFolder'}."/".$configHash->{'ourRiboSeqs'}) 
	   and (-e $configHash->{'GFFAndOutputFolder'}."/".$configHash->{'R_seqs_stats_script'})
	   and (-e $configHash->{'GFFAndOutputFolder'}."/".$configHash->{'R_plot_script'})
        and (-d $configHash->{'DBCreationFolder'}) 
          and (-d $configHash->{'ProgExecFolder'}) 
            and (-d $configHash->{'GFFAndOutputFolder'}) 
              and (-d $configHash->{'usefulFolder'}) 
                and (-d $configHash->{'guideFolder'}) 
                  and (-d $configHash->{'configFolder'}) ) ){      
    annoDie( "ERROR: This Annocript version is corrupted. Please download and re-install Annocript. Exiting..\n");
  }    
  
  #Checking folders and files names
  #Names of the output of the programs
  if ( $configHash->{'rpstblastnOut'} !~ /^[A-za-z0-9\_\-]*$/i  ){
		annoDie( "The name".$configHash->{'rpstblastnOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'blastxSPOut'} !~ /^[A-za-z0-9\_\-]*$/i  ){ 
  
		annoDie( "The name".$configHash->{'blastxSPOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'blastxTROut'} !~ /^[A-za-z0-9\_\-]*$/i  ){
		annoDie( "The name".$configHash->{'blastxTROut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'blastxUnirefOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'blastxUnirefOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'blastnOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'blastnOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'portraitOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'portraitOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'dna2pepOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'dna2pepOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'outCount'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'outCount'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  #Names of the stats files
  if ( $configHash->{'bpStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'bpStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'mfStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'mfStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'ccStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'ccStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'cdStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'cdStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
    if ( $configHash->{'pathwaysL1File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'mfStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'pathwaysL2File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'ccStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'pathwaysL3File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		annoDie( "The name".$configHash->{'cdStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  #Filenames with a dot
  if ( $configHash->{'outFileName'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'outFileName'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'outFiltered'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'outFiltered'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'ORFFastaFileName'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'ORFFastaFileName'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'NCOutFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'NCOutFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'gcContentFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'gcContentFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'closerOSFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'closerOSFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  if ( $configHash->{'ATGCPercFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'ATGCPercFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  if ( $configHash->{'outHashFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'outHashFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
 
  if ( $configHash->{'parsingFilePath'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'parsingFilePath'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  if ( $configHash->{'filtMappingFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'filtMappingFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  if ( $configHash->{'R_plot_script'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'R_plot_script'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}  
  
  if ( $configHash->{'R_seqs_stats_script'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'R_seqs_stats_script'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
  }  
	
  if ( $configHash->{'R_log_file'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'R_log_file'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  if ( $configHash->{'gffDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'gffDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	} 
  
    if ( $configHash->{'ourRiboSeqs'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'ourRiboSeqs'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}  
  
  if ( $configHash->{'rfamDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'rfamDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}
  
  if ( $configHash->{'gffDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'gffDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n");
	}   
  
  #Filenames with HTML extension
  if ( ($configHash->{'htmlHome'} !~ /^[A-za-z0-9\_\-\.]*$/i ) and ($configHash->{'htmlHome'} !~ /.html/ )){
	annoDie( "The name".$configHash->{'htmlHome'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n");
  }   
  if ( ($configHash->{'htmlAnnotStats'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'htmlAnnotStats'} !~ /.html/ ) ){
	annoDie( "The name".$configHash->{'htmlAnnotStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n");
  }   
  if ( ($configHash->{'htmlSeqStats'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'htmlSeqStats'} !~ /.html/ ) ){
	annoDie( "The name".$configHash->{'htmlSeqStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n");
  }
	
 #Logo name check
  if ( ($configHash->{'AnnocriptLogo'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'AnnocriptLogo'} !~ /.jpeg|.png/ ) ){
	annoDie( "What will you use as logo?! ".$configHash->{'htmlSeqStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the png or jpeg format. Program will exit...\n");
  }
		  
  #Checking other variables 
  checkVariable($configHash->{'keepTempFiles'},'keepTempFiles',"Temporary files will be stored\n");
 
  if (!correct_type($configHash->{'maxLengthDescs'},"positiveint")){
    annoDie( $configHash->{'maxLengthDescs'}." is not a positive integer. Please change it and re-execute Annocript. Exiting...\n");
  } 
  
  if (!correct_type($configHash->{'bigMemory'},"positiveint")){
    annoDie( $configHash->{'bigMemory'}." is not a positive integer. Please change it and re-execute Annocript. Exiting...\n");
  }
}



=head2 check_query_file

 Title   : check_query_file
 Usage   : check_query_file(  - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function:   check if the query FASTA file is well written and is in the working folder. Behaviour of the function is different
          if user is calling a saved session. The file is moved soon in the session folder. 
          
 Returns : nothing

=cut
sub check_query_file {

#CONTROLS IF THE QUERY EXISTS AND THE NAME IS .fa or .fasta
  my $fastaSeqs = $configHash->{'fastaSeqs'};
  
  #The name has to be without dots inside. Only one dot can separate the extension from the name
  my $dots=0;
  my @line = split (//,$fastaSeqs);
  foreach my $char(@line){
      if ($char eq '.'){
          $dots++;
      }
  }
  if ($dots>1){
      annoDie( "Please change the name of your file removing all internal dots (not the one used for the extension). Annocript will exit...\n");
  }
  
  #annoPrint ("Fasta name read: ".$fastaSeqs."\n";
  my @name = split(/\./,$fastaSeqs);
  #annoPrint ("Fasta name read: ".$name[0]."\n";
  #Here is created the name to use for all the files of output by extracting it from the fasta file name
  #It cannot be longer than 50 chars
  if (@name > 1){
    #$configHash->{'name'} = extract_name($configHash->{'sessionFolder'},0);
    $configHash->{'name'} = $name[0]; 
    #Checks for permitted characters
    if ( (length($configHash->{'name'}) > 50) # or !($configHash->{'fastaSeqs'} =~ /(\.fasta?)|(\.fa?)/) 
            or ($configHash->{'fastaSeqs'} !~ /^[A-za-z0-9\_\-]+\.(fa|fasta)$/i)){
            annoDie( $configHash->{'fastaSeqs'}." is not a correct FASTA file name. Allowed characters [A-Z,a-z,0-9,_,-]. Allowed extensions [.fasta,.fa]. Max length: 50 chars\n");
    }else{
      
      #If session is a new one
      #Fasta file must stay only in the working folder and it will be shifted. If it is not there,then only db creation can be done or die
      if ( $newSession == 1 ){
        my $fastaPath = $workingFolder.'/'.$configHash->{'fastaSeqs'};
        annoPrint ("Path to control for fasta: $fastaPath\n");
        #The Sequences file can also stay already in the session folder created by the user in a previous execution
        if (-e $fastaPath) {
          annoPrint ("Sequences found in your working folder. Copying in session folder...\n");
          #move($fastaPath,$sessionFolder);
          copy($fastaPath,$sessionFolder) or annoDie("Cannot copy $fastaPath in $sessionFolder");
          $configHash->{'fastaSeqsPresent'} = 'YES';
        }else{
          #This else happens when sequences file is not there. Then, if the user chose to execute analysis or to print output
          # Annocript has to die.
          if ( ($configHash->{'doExecutePrograms'} eq 'YES') or ($configHash->{'doBuildOutput'} eq 'YES')
            or ($configHash->{'extractStatistics'} eq 'YES')){
              annoDie( "No sequences file found in $fastaPath. You'll cannot execute analysis...\n");
            }
          annoPrint ("No sequences file found in  $fastaPath. But you are only creating a DB...\n");
          $configHash->{'fastaSeqsPresent'} = 'NO';
        }
      }else{#If session is loaded then the fasta must stay only in the session folder. No other fasta can be shifted there. Die otherwise.
        if (-e $sessionFolder.'/'.$configHash->{'fastaSeqs'}) {
          annoPrint ("Sequences found in your session folder...\n");
          $configHash->{'fastaSeqsPresent'} = 'YES';
        }else{
          annoDie("Unable to find ".$configHash->{'fastaSeqs'}." in $sessionFolder. Annocript will go and you can figure out why.".
                    " \nThere are 2 possible reasons: \n\t 1. you written a wrong fasta name in the configuration file.".
                    "\n\t 2. the name in the config file is right but not the folder you chose\n Bye!\n ");
        }
      }
      
      if ($configHash->{'fastaSeqsPresent'} eq 'YES'){ 
        #Here we control what kind of alphabet has the sequence and store the value in a variable
        $configHash->{'sequencesAlphabet'} = detect_fasta_type($sessionFolder.'/'.$fastaSeqs);
        annoPrint ("Annocript detected that your sequences are ");
        if ( ($configHash->{'sequencesAlphabet'} eq 'dna') or ($configHash->{'sequencesAlphabet'} eq 'rna') ){
            annoPrint ("transcripts\n");
        }else{
            annoPrint ($configHash->{'sequencesAlphabet'}."\n");
         }
        #Here we count and store the number of sequences in the fasta file
        my $seqsPath = $sessionFolder.'/'.$fastaSeqs;
        my $totSequences = `fgrep -c '>' $seqsPath`; #Extract number of sequences from the query file
        annoPrint("Number of sequences: $totSequences");
        $configHash->{'totSequences'} = $totSequences;
      }
    }
  }else {
    annoDie( "The name of file with transcripts needs to have a .fa or .fasta extension. Exiting...\n");
  }
}



=head2 checkPROGRAMS_EXEC

 Title   : checkPROGRAMS_EXEC
 Usage   : checkPROGRAMS_EXEC( - configHash = the piece of the config hash for PROGRAMS_EXEC
							 - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function:  Checks the config hash parameters regarding the execution of programs: how they are written, if they are present and so on..
 Returns : nothing

=cut
sub checkPROGRAMS_EXEC {
  my $localFolder = shift;
  
  #annoPrint (Dumper\$configHash;
	#CONTROL THE PATH TO DATABASES  
	if ( $configHash->{'rfamDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'rfamDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n");
	}   
	if (  $configHash->{'swissProtDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'swissProtDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n");
	}
	if ( $configHash->{'tremblDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		annoDie( "The name".$configHash->{'tremblDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n");
	}
             


  
  #Check if the links wrote for the programs are valid
  check_programs_presence();     
  
  #Check the blast input parameters using the final chars near the parameter
  #X is used for BLASTX and BLASTP, N for BLASTN and TBLASTN, RPS for RPSTBLASTN and RPSBLAST
  my @blastIds = ('X','N','RPS'); 
          
  foreach my $blastId (@blastIds){
    my $arg = 'word_size'.$blastId;
			if (defined($configHash->{$arg})){  
				 if (!correct_type($configHash->{$arg},"positiveint")){
          annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
          }  
      }
      $arg = 'evalue'.$blastId;
      if (defined($configHash->{$arg})){  
        if (!correct_type($configHash->{$arg},"real")){
          annoDie( $configHash->{$arg}." is not a real. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}  
			}
			$arg = 'num_descriptions'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}  
			}
			$arg = 'num_alignments'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{ $arg},"positiveint")){
					annoDie( $configHash->{ $arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}  
			}
			$arg = 'max_target_seqs'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}  
			}
      $arg = 'num_threads'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}else{
          my $num_processors = num_processors();
       
          if ($configHash->{$arg} > $num_processors ) {
              annoDie ("ERROR: Please change the paramater $arg in config file. ".
              " Keep in mind that your processors are: $num_processors.\n");
          }
        }  
			}
      $arg = 'threshold'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}  
			}
      $arg = 'outfmt'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					annoDie( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n");
				}else{
            
        }  
			}
		}
    
    #Number of processes to use for the parallel execution 
    #This checks if your choice is good depending on the number of processes used      
		if (defined($configHash->{'threads4Parallel'})){  
			if (!correct_type($configHash->{'threads4Parallel'},"positiveint")){
				annoDie( $configHash->{'threads4Parallel'}." is not a positive integer. Please change threads4Parallel and re-execute Annocript. Exiting...\n");
			}else { 
        my $num_processors = num_processors();
       
        if ($configHash->{'threads4Parallel'} > $num_processors ) {
            annoDie ("ERROR: Please change the paramater threads4Parallel in config file. ".
            " Keep in mind that your processors are: $num_processors.\n");
        }      
      }  
		}
    	 
		# Chekcing correctness of PORTRAIT's parameters
		if (defined($configHash->{'NCSeqLength'})){  
			if (!correct_type($configHash->{'NCSeqLength'},"positiveint")){
				annoDie( $configHash->{'NCSeqLength'}." is not a positive integer. Please change NCSeqLength and re-execute Annocript. Exiting...\n");
			}  
		}
		if (defined($configHash->{'NCORFLength'})){  
			if (!correct_type($configHash->{'NCORFLength'},"positiveint")){
				annoDie( $configHash->{'NCORFLength'}." is not a positive integer. Please change NCORFLength and re-execute Annocript. Exiting...\n");
			}  
		}
		if (defined($configHash->{'NCThresh'})){  
			if (!correct_type($configHash->{'NCThresh'},"real")){
				 annoDie( $configHash->{'NCThresh'}." is not a real. Please change NCThresh and re-execute Annocript. Exiting...\n");
			}  
		}
		
    if (defined($configHash->{'max_port_seq_length'})){  
			if (!correct_type($configHash->{'max_port_seq_length'},"positiveint")){
				annoDie( $configHash->{'max_port_seq_length'}." is not a positive integer. Please change max_port_seq_length and re-execute Annocript. Exiting...\n");
			}  
		}
    if (defined($configHash->{'min_port_seq_length'})){  
			if (!correct_type($configHash->{'min_port_seq_length'},"positiveint")){
				annoDie( $configHash->{'min_port_seq_length'}." is not a positive integer. Please change min_port_seq_length and re-execute Annocript. Exiting...\n");
			}  
		}
    #Checking parameters to use for dna2pep
    if (defined($configHash->{'d2pMode'})){  
			if ( ($configHash->{'d2pMode'} ne "strict") and ($configHash->{'d2pMode'} ne "any") and ($configHash->{'d2pMode'} ne "none") ){
				 annoDie(" d2pMode should be one among: strict, any, none. And it is: ". $configHash->{'d2pMode'}.". Please change it and re-execute Annocript. Exiting...\n");
			}  
		}    
      
    
    #Annocript Checks if there is a file of counts and if it is good
    #check_counts_file();
    #IF ONE DAY WE'LL DECIDE TO USE THE COUNTS WE CAN START FROM HERE... TILL THAT...
    $configHash->{'countOutPresent'} = 'NO';
}


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args( );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

  #  Parse options
  GetOptions(
           "help",         \$HELP,
           "version",      \$VERSION
          );
  #Print a little help  
  if ($HELP){
    pod2usage(1);
    exit;
  }
  #Print version
  if ( $VERSION ){
    print "Annocript.pl release $version \n";
    exit;
  }
  
  ## If no arguments were given, then allow STDIN to be used only
  ## if it's not connected to a terminal (otherwise annoPrint (usage)
  pod2usage("$0: Config file is missing...")  if ( (@ARGV == 0) && (-t STDIN));
   
  #Exit if no config file is inputed  
  if (! (-e $ARGV[0]) ) {
    exit "Config file does not exists. Please write a correct path.\n";
  }
}
