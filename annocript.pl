#!/usr/bin/perl

#Annocript - A complete tool for transcriptomes annotation
#Copyright (C) <2015>  <Francesco Musacchia>

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

Annocript is a tool for the annotation of transcriptomes.

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

 $Id: Annocript, 2015/14/19 11:00:12 Exp $
=cut

use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Cwd;#To change work directory
 
#USE THE GENERAL UTILITY MODULES
use File::Copy;#To manage file
use FindBin;#To search libraries
use lib $FindBin::Bin;#To search libraries
use Term::ReadKey;
use Getopt::Long;#To control the input
use Pod::Usage;#Used to write an usage
use LWP::Simple;#Functions for FTP management
use File::Path qw(make_path remove_tree);#To remove directories
use IO::Handle;#To immediately print with autoflush 
 
#Annocript libraries
use LIB::db_creator qw( execute_db_creator );
use LIB::execute_programs qw( execute_programs  );
use LIB::programs_out_parser qw(execute_programs_out_parser);

use LIB::programs_management qw(log_and_exit print_and_log checkVariable initialize_folders
															checkConfigVariables configFile2Hash);

use LIB::annocript_management qw(check_annocript_integrity check_input_fasta_file
									valid_analysis check_uniref_presence check_programs_presence
									uniprot_db_selection detect_gff_db check_db_tables_version
									getDBVersion clean_cd_name nice_printing indexed_db_present);
									
use LIB::html_management qw(testFTP	ftp_fetch_file check_FTP_diff_sources checkLink
														check_url);

use LIB::mysql_db_management qw(check_DB_user_permissions db_present);

#Using a library to manage files
use LIB::files_management qw(  extract_name check_presence);
									
#Using a library for standard utilities								
use LIB::std_lib qw(num_processors ram_memory correct_type save_hash);
		
use LIB::utils qw(execute_utils);

#EXTRACT PARAMETERS FROM THE CONFIG FILE
my $configHash;#This hash holds up all the life of Annocript
my $program_name = "annocript.pl";
my $version = '2.0';
my $dbTablesVersion = '0.3';
my $user_config = "";
my $dataFolder = "data";
my $logFolder = "log";
my $lib_folder = "LIB";
my $analysesFolder = "analyses";
my $outFolder ="output";
my $gffFolder = "gff";
my $statsFolder = "stats";
my $tempTableFolder = 'tables';#Folder to keep the tables while they are being generated

my $usefulFolder = 'USEFUL';#A folder with useful scripts and data
my $workingFolder = '';#Folder where Annocript has to work
my $foldersFile = "folders.txt";
my $log_file = "";
my $timesFile = "";

#Variables present also in utils.pl (IF YOU CHANGE IT HERE CHANGE ALSO IN UTILS.PL!!)
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



#NEW VARIABLES
my $config_folder = "CONFIGURATION";
my $program_config = "program_config.txt";#The path to the configuration file for the program
my $license_file = "COPYRIGHT";
my $license = "";
my $variables_file = "variables.txt";#The file with all the variables that must be present in the config files
my $analysisFolder = "";
my $analysis = "";
my $session = "";
my $exec_utils = 0;
my $utils = "";
my $program_folder = "";
my $working_folder = "";
my $uniref_db_ver = "";


my $force_run = "";

my $newSession = 0;

#Database selection
my $db_version = "";
my $db_type = "";

my $niceClosure = "\n=================================================\n";
my $configClosure = "\n*************************************************\n";

my $program_title = "***********************************************************************\n".
		"This program comes with ABSOLUTELY NO WARRANTY; for details http://opensource.org/licenses/GPL-3.0.
This is free software, and you are welcome to redistribute it
under certain conditions; http://opensource.org/licenses/GPL-3.0 for details.\n".
"***********************************************************************\n".
"  /&&&&&&                                              /&&           /&&          \n".
" /&&__  &&                                            |__/          | &&          \n".
"| &&  | &&/&&&&&&& /&&&&&&&   /&&&&&&  /&&&&&&& /&&&&&& /&& /&&&&&& /&&&&&&       \n".
"| &&&&&&&| &&__  &| &&__  && /&&__  &&/&&_____//&&__  &| &&/&&__  &|_  &&_/       \n".
"| &&__  &| &&  | &| &&  | &|  &&  | &| &&     | &&  |__| &| &&  | && | &&         \n".
"| &&  | &| &&  | &| &&  | &|  &&  | &| &&     | &&     | &| &&  | && | && /&&     \n".
"| &&  | &| &&  | &| &&  | &|   &&&&&&|  &&&&&&| &&     | &| &&&&&&&/ |  &&&&/    \n".
"|__/  |__|__/  |__|__/  |__/ |______/ |_______|__/     |__| &&____/   |___/      \n". 
"                                                          | &&\n".                        
"                                                          | &&\n".                                           
"                                                          |__/\n".  
" \nGiven a fasta file with a transcriptome, Annocript will annotate your sequences and".
" it will separate putative coding and long non-coding RNA sequences.\n".
" Reference: Musacchia F, Basu S, Petrosino G, Salvemini M, Sanges R. \n".
" Annocript: a flexible pipeline for the annotation of transcriptomes which can also identify ".
" putative long non-coding RNA.".
" Bioinformatics (2015) doi:10.1093/bioinformatics/btv106".
"\n***********************************************************************\n\n\n";
											
											
											
######################################################### MAIN ########################################################


#LOG FILE CREATION
#Here we take the current time to get a log name different at each computation  
my $time = scalar(localtime);
$time =~ s/ /_/g;
$time =~ s/:/-/g;

STDOUT->autoflush(1);#This makes STDOUT hot in the sense that everything will be print immediately


#Annocript permits both the usage of some utilities and the run of an analysis pipeline
#we set here the log file name depending on what is doing
if ( $utils ){
	$log_file = "utils_$time.log";
}else{
	$log_file = "annocript_exec_$time.log";
	$timesFile = "annocript_times_$time.log";
}

#open the log file
my $logFHandle;
open($logFHandle, ">$log_file") or die "ERROR [$!]: Cannot open $log_file! Check permissions.\n";


#  Parse command line arguments.  
parse_command_line_args();




#If there is no utility to do, go ahed and check the input configuration file
if ( -e $user_config and $user_config ne ''){
	print_and_log("\n\nPRELIMINARY CONFIGURATION...\n",$log_file);#DEBUGCODE
	
	#Sets useful directories and stores the content of the configuration files
	set_global_var_and_confighash($user_config);

	#Checks folder and files if they are well written
	print_and_log("Checking Annocript integrity...\n",$log_file);#DEBUGCODE
	check_annocript_integrity($configHash,$log_file);
		
	#Print license
	if ( $license ){
		print_file($program_folder."/".$license_file);
		log_and_exit("Exiting..\n",$log_file);
	}  
		
	#Utilities execution
	if ( $utils ){
		print_and_log("Executing utilities...\n",$log_file);#DEBUGCODE
		execute_utils($program_config);		
	}else { 
	
	print_and_log("Checking input parameters...\n",$log_file);#DEBUGCODE	
	#Detect the analysis to use
	check_input_parameters($analysis);	
	
	#Now that we have all the parameters and the analysis name
	#I check all the parameters and configuration 		
	configure_analysis();

	#Saves the hash in a file to be used later
	my $hashName = "configHash";
	my $configHashPath = $dataFolder."/".$hashName;
	save_hash($configHashPath,$configHash, $hashName);  
	my $newLogFile =  $logFolder."/".$log_file;#Log of Annocript goes in the log folder

	print_and_log( "Your configuration has been saved! A log file will be created in $logFolder/$log_file. \nThe overall computation usually"
		  ." depends on the number of sequences and the speed of the machine you are using. \nYou may want to check Annocript step-by-step "
		  ."with the following command: more ".$newLogFile.".\n"
		  ."\nYou will find all the results in the folder you chose here. Please take a look at\n"
		  ." https://github.com/frankMusacchia/Annocript/blob/master/GUIDE/OUTPUT.md to understand the output from Annocript and the"
		  ." organization of folders."
		  ."\n\nPlease let us know if any problem occurs during the computation or you don't get from Annocript the expected result.\n"
		  ."You may want to use the forum at: https://groups.google.com/forum/#!forum/annocript.\n",$log_file);#DEBUGCODE	

	nice_printing("STARTING ANNOCRIPT IN BACKGROUND! BYE!",$log_file);

	close ( $logFHandle );


	#print "moving the file now...";
	move($workingFolder.'/'.$log_file, $newLogFile) or log_and_exit(" Annocript will not move $workingFolder/$log_file to $newLogFile\n");
	$log_file = $newLogFile;
	$timesFile = $logFolder."/".$timesFile;#File with times goes inside the log folder


	#We also remove all the other remaining log file in the folder created by user's interruptions
	system("rm -f $workingFolder/annocript_exec_*"); 
		

	#Copying the user config from working to analysis folder
	(system("cp $user_config $analysisFolder") ) == 0 
	  or log_and_exit("ERROR: Unable to copy $user_config in $analysisFolder. Please check permissions...\n",$log_file);              	
	#Make the user config complete path
	$user_config = $analysisFolder."/".$user_config;
	#Go to program folder directory
	chdir $program_folder;
	
	#This perl script will be run in background
	my $runAnnocriptCmd = "nohup perl ".$configHash->{'ProgExecFolder'} ."/annocript_executor.pl $user_config $program_config $timesFile $configHashPath >> $log_file &";
	( system($runAnnocriptCmd) ) == 0 
	  or die("Unable to start annocript_executor!\n");
	}
}else{
	log_and_exit("Cannot open configuration file. Please give a correct parameter for user_config! Exiting...\n",$log_file);
}


#########################################NEEDED ROUTINES

=head2 configure_analysis

 Title   : configure_analysis
 Usage   : configure_analysis( - configFilePath = path of the config file
                               );

 Function:  this subroutine completes the configuration of annocript:
              - checks the variables with YES, NO values
              - creates needed folders
              - check variables related to the three modules (DB_CREATOR, PROGRAMS_EXEC, PROGRAMS_OUT_PARSER)
              - build some important variables in the configuration hash
 Returns : nothing

=cut
sub configure_analysis{
   
	#Here I put a control of the variables in the configuration file that can be YES or NO
	print_and_log ("\nAnnocript is going to run:\n",$log_file);
	checkVariable($configHash->{'doDbCreation'},'doDbCreation', "\t-Database Creation.\n");
	checkVariable($configHash->{'doExecutePrograms'},'doExecutePrograms', "\t-Execution of programs:\n");

	if ($configHash->{'doExecutePrograms'} eq 'YES'){
	checkVariable($configHash->{'doBlastxSP'},'doBlastxSP',"\t\tBLASTX against SwissProt.\n");
	checkVariable($configHash->{'doBlastxTRorUf'},'doBlastxTRorUf',"\t\tBLASTX against TrEMBL/UniRef.\n");
	checkVariable($configHash->{'doRpstblastn'},'doRpstblastn',"\t\tRPSBLAST against conserved domains database.\n");
	checkVariable($configHash->{'doBlastn'},'doBlastn',"\t\tBLASTN against ncRNA.\n");
	checkVariable($configHash->{'doLNCPrediction'},'doLNCPrediction',"\t\tlncRNA prediction. \n");
	checkVariable($configHash->{'doDna2Pep'},'doDna2Pep',"\t\tLongest ORF search.\n");
	}
	checkVariable($configHash->{'doBuildOutput'},'doBuildOutput', "\tCreation of the final tabular output.\n");
	checkVariable($configHash->{'extractStatistics'},'extractStatistics',"\tGeneration of statistics.\n\n");

	checkVariable($configHash->{'useGFFDB'},'useGFFDB',"\tGFF3 database will be generated.\n");
	checkVariable($configHash->{'printGFFOutput'},'printGFFOutput',"\tGFF3 files will be printed.\n");


	#Checks everything regarding the query file. If session is new, the file is put in the session folder
	$configHash->{'name'} = check_input_fasta_file($configHash->{'fastaSeqs'},\$configHash,$newSession,$workingFolder,$analysisFolder);

	$dataFolder = $analysisFolder."/".$configHash-> {'dataFolder'};
	#Check if directory DATA exists, otherwise it creates it
	unless(-d $dataFolder){
		print_and_log ("Creating folder $dataFolder...\n",$log_file);
		mkdir $dataFolder or log_and_exit("ERROR: can't create folder $dataFolder. Check permissions. \n",$log_file);
	}

	$logFolder = $analysisFolder."/".$configHash-> {'logFolder'};
	#Check if directory LOG exists, otherwise it creates it
	unless(-d $logFolder){
		print_and_log ("Creating folder $logFolder...\n",$log_file);
		mkdir $logFolder or log_and_exit("ERROR: can't create folder $logFolder. Check permissions. \n",$log_file);
	}

	$outFolder = $analysisFolder."/".$configHash->{'outFolder'};
	#Check if directory OUTPUT exists, otherwise it creates it
	unless(-d $outFolder){
	print_and_log ("Creating folder $outFolder....\n",$log_file);
	mkdir $outFolder or log_and_exit("ERROR: can't create folder $outFolder. Check permissions. \n",$log_file);
	}

	$statsFolder = $analysisFolder."/".$configHash->{'statsFolder'};#Statistics folder
	#Creation of a folder for the statistics
	unless(-d $statsFolder){
		print_and_log ("Creating folder $statsFolder...\n",$log_file);
		mkdir $statsFolder or log_and_exit("ERROR: can't create folder $statsFolder. Check permissions. \n",$log_file);
	}#else{print_and_log ("Folder $statsFolder already exists. Continuing...\n ",$log_file);}

	$gffFolder = $analysisFolder."/".$configHash->{'gffFolder'};
	#Check if directory for GFF files exists, otherwise it creates it
	unless(-d $gffFolder){
		print_and_log ("Creating folder ".$gffFolder."...\n",$log_file);
		mkdir $gffFolder or log_and_exit( "ERROR : can't create folder $gffFolder. Check permissions. \n",$log_file);
	} 

	#Check if directory for BLASTX output exists, otherwise it creates it
	unless(-d $outFolder."/blastx"){
		print_and_log ("Creating folder ".$outFolder."/blastx...\n",$log_file);
		mkdir $outFolder."/blastx" or log_and_exit( "ERROR : can't create folder $outFolder /blastx. Check permissions. \n",$log_file);
	} 
	#Check if directory for RPSTBLASTN output exists, otherwise it creates it
	unless(-d $outFolder."/rpstblastn"){
		print_and_log ("Creating folder ".$outFolder."/rpstblastn...\n",$log_file);
		mkdir $outFolder."/rpstblastn" or log_and_exit( "ERROR : can't create folder $outFolder /rpstblastn . Check permissions. \n",$log_file);
	}  
	#Check if directory for BLASTN output exists, otherwise it creates it
	unless(-d $outFolder."/blastn"){
		print_and_log ("Creating folder ".$outFolder."/blastn...\n",$log_file);
		mkdir $outFolder."/blastn" or log_and_exit( "ERROR : can't create folder $outFolder /blastn. Check permissions. \n",$log_file);
	}

	  
	#print_and_log (Dumper\$configHash,$log_file);    #DEBUGCODE
	
	$configHash->{'force_run'} = $force_run;
	#The version of Annocript
	$configHash->{'AnnocriptVer'} = $version;
	$configHash->{'dbTablesVersion'} = $dbTablesVersion;
	$configHash->{'usefulFolder'} = $usefulFolder;
	$configHash->{'faSomeRecords'} = $program_folder."/$usefulFolder/".$configHash->{'faSomeRecords'};
	

	#print_and_log ("\n>>> Database settings (it will need a while):\n",$log_file);
	checkDB_CREATION();
	#print_and_log ("\n>>> Programs settings:\n",$log_file);
	checkPROGRAMS_EXEC();
	#print_and_log ("\n>>> Setting final statistics and outputs\n",$log_file);
	checkPROGRAMS_OUT_PARSER();
	print_and_log ("...DONE!\n",$log_file);

	my $dbDataFolder = $configHash->{'dbDataFolder'};

	#Adding some variables for blast out types. This is done 
	#because if a day we want to add more types, nothing is to change in gff3_manager.pm
	$configHash->{'blastxSPOutType'} = $configHash->{'blastOutType'};
	$configHash->{'blastxTROutType'} = $configHash->{'blastOutType'};
	$configHash->{'blastxUnirefOutType'} = $configHash->{'blastOutType'};
	$configHash->{'rpstblastnOutType'} = $configHash->{'blastOutType'};
	$configHash->{'blastnOutType'} = $configHash->{'blastOutType'};

	#Writing the complete path (with datafolder) of the databases    
	$configHash->{'ncDB'} = "$dbDataFolder/".$configHash->{'ncDB'};

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

	#Check dependencies
	check_dependencies();
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
    if( ((scalar <$outFolder/rpstblastn/*>) or (scalar <$outFolder/blastx/*>) or (scalar <$outFolder/blastn/*>)) and ($force_run eq '')) {
     log_and_exit("\n WARNING: Some output files already exist! If you will run a program (BLAST, lncRNAPrediction, DNA2PEP) the corresponding ".
                                "output must be overwritten. Please use the --force_run parameter then.... Exiting\n",$log_file);
    }
  #If programs will not be run and the user wants to print the final output table then a check is needed if outputs of programs
  #are present
  }elsif ( ($configHash->{'doExecutePrograms'} ne 'YES') and ($configHash->{'doBuildOutput'} eq 'YES') ){
    if( !(scalar <$outFolder/rpstblastn/*>) and !(scalar <$outFolder/blastx/*>) and !(scalar <$outFolder/blastn/*>)) {
      log_and_exit("Output folders are empty. The output from Annocript cannot be obtained without that.".
      "Please chose doExecutePrograms = YES in the configuration file. Exiting...\n",$log_file);
    }
  #If only db creation is active we only check the presence of the makeblastdb program
  }elsif ( $configHash->{'doDbCreation'} eq 'YES'){
    if (!check_presence($configHash->{'makeblastdbPath'})){
        log_and_exit("You miss makeblastdb. Please check the links in the configuration file or re-install the BLAST program. Exiting...\n",$log_file);
    }
  } 
  
  #If nothing is there in the folders we cannot start the output and stats creation
  if ( ($configHash->{'extractStatistics'} eq 'YES') and (($configHash->{'doBuildOutput'} ne 'YES')) ){
            if (!check_presence($outFolder.'/'.$configHash->{'outFiltered'})){
                log_and_exit("You cannot start the statistics without the filtered ouput of Annocript. Please restart with doBuildOutput = YES\n",$log_file);
            }
  }
  
}


###############################Modules Execution checks

=head2 checkDB_CREATION

 Title   : checkDB_CREATION
 Usage   : checkDB_CREATION( - configHash = the piece of the config hash for DB_CREATION
							 - local folder = the folder where actually is the DB_CREATION script
                               );

 Function:   Checks the config hash parameters regardin the database creation and user. How they are written, if they are present and so on..
 Returns : nothing

=cut
sub checkDB_CREATION {

  my $wrongLinks = 0;  
 
  ##########################CHECKING THE LINKS   
  if ($configHash->{'doDbCreation'} eq 'YES'){
    
    #Check permissions. If the user does not have, he/she cannot create database or perform analyses
    if (!check_DB_user_permissions($configHash->{'platform'},$configHash->{'host'},
              $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'})) {
        log_and_exit("ERROR: You do not have enough privileges to create the DB. Please ask your adiministrator. Annocript will close...\n",$log_file);
    }
      # print_and_log (Dumper\$configHash,$log_file);
      
      #A check on email address 
      if ( !($configHash->{'uniprotWebPass'} =~ /@/) or !($configHash->{'uniprotWebPass'} =~ /\.[a-zA-Z]+?/) ){
        log_and_exit( "ERROR : ".$configHash->{'uniprotWebPass'}." is not a correct email address. Please check it in config file \n",$log_file);
      }
      
      #Triple controls on the links. We need it three times because it can fail.
      print_and_log ("Checking the links for the databases to download...(it may take a while)...\n",$log_file);
			
			if (!check_url($configHash->{'swissprotDBLink'}) ) {
         my $new_link = check_FTP_diff_sources($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
          if ( $new_link ne ''){
              $configHash->{'swissprotDBLink'} = $new_link;
         }else{
            print_and_log ("Please check if [".$configHash->{'swissprotDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
            $wrongLinks++;
           }
         }

        if (!check_url($configHash->{'tremblDBLink'}) ) {
          my $new_link = check_FTP_diff_sources($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
            if ( $new_link ne ''){
                $configHash->{'tremblDBLink'} = $new_link;
            }else{
              print_and_log ("Please check if [".$configHash->{'tremblDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
              $wrongLinks++;
             }
         }
 
         if (!check_url($configHash->{'uniprotVerLink'}) ) {
          my $new_link = check_FTP_diff_sources($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
          if ( $new_link ne ''){
              $configHash->{'uniprotVerLink'} = $new_link;
          }else{
             print_and_log ("Please check if [".$configHash->{'uniprotVerLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
            $wrongLinks++;
           }
        }
         if (!check_url($configHash->{'unirefDBLink'}) ) {
          my $new_link = check_FTP_diff_sources($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
          if ( $new_link ne ''){
              $configHash->{'unirefDBLink'} = $new_link;
          }else{
             print_and_log ("Please check if [".$configHash->{'unirefDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
            $wrongLinks++;
           }
         }       
 
         if (!check_url($configHash->{'unirefVerLink'}) ) {
          my $new_link = check_FTP_diff_sources($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
          if ( $new_link ne ''){
              $configHash->{'unirefVerLink'} = $new_link;
          }else{
            print_and_log ("Please check if [".$configHash->{'unirefVerLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
            $wrongLinks++;
           }
         }
 
         if (!check_url($configHash->{'GODBLink'}) ) {
          my $new_link = check_FTP_diff_sources($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
          if ( $new_link ne ''){
              $configHash->{'GODBLink'} = $new_link;
          }else{
            print_and_log ("Please check if [".$configHash->{'GODBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
            $wrongLinks++;
           }
         }
         
         if (!check_url($configHash->{'enzymeDBLink'}) ) {
          print_and_log ("Please check if [".$configHash->{'enzymeDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file);#
          $wrongLinks++;
         }
         if (!check_url($configHash->{'cdDBLink'})) {
          annoPrint ("Please check if [".$configHash->{'cdDBLink'}."] is a correct URL. Annocript will continue...\n"); #
          $wrongLinks++;
         }        
         if (!check_url($configHash->{'cdTableLink'})) {
           print_and_log ("Please check if [".$configHash->{'cdDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
          $wrongLinks++;
         }
         if (!check_url($configHash->{'ncRNADBLink'})) {
          print_and_log ("Please check if [".$configHash->{'ncRNADBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
          $wrongLinks++;
         }        

         if (!check_url($configHash->{'pathwaysTableLink'})) {
          print_and_log ("Please check if [".$configHash->{'pathwaysTableLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
          $wrongLinks++;
         }        
                                                              
      ##This IFs contain the check of the links to the databases. Two types of check are present for links to the Uniprot data.
       ##LWP::Simple::Head and NET::Ftp
      ##print "Swiss-prot DB link...";
        #if (!head($configHash->{'swissprotDBLink'}) ) {
          #if (testFTP($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #my $new_link = check_FTP_diff_sources($configHash->{'swissprotDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
              #if ( $new_link ne ''){
                  #$configHash->{'swissprotDBLink'} = $new_link;
              #}else{
                #print_and_log ("Please check if [".$configHash->{'swissprotDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
                #$wrongLinks++;
              #}
            #}
          #}
        #}
       ##print "Trembl DB link...";
        #if (!head($configHash->{'tremblDBLink'}) ) {
          #if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #if (testFTP($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #my $new_link = check_FTP_diff_sources($configHash->{'tremblDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                #if ( $new_link ne ''){
                    #$configHash->{'tremblDBLink'} = $new_link;
                #}else{
                  #print_and_log ("Please check if [".$configHash->{'tremblDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
                  #$wrongLinks++;
                #}
              #}
            #}
          #}
        #}
        ##print "Uniprot DB version link...";
        #if (!head($configHash->{'uniprotVerLink'}) ) {
          #if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #if (testFTP($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                #my $new_link = check_FTP_diff_sources($configHash->{'uniprotVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                #if ( $new_link ne ''){
                    #$configHash->{'uniprotVerLink'} = $new_link;
                #}else{
                  #print_and_log ("Please check if [".$configHash->{'uniprotVerLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
                  #$wrongLinks++;
                #}
              #}
            #}
          #}
        #}  
        ##print "UniRef DB link...";
        #if (!head($configHash->{'unirefDBLink'}) ) {
          #if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #if (testFTP($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                #my $new_link = check_FTP_diff_sources($configHash->{'unirefDBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                #if ( $new_link ne ''){
                    #$configHash->{'unirefDBLink'} = $new_link;
                #}else{
                  #print_and_log ("Please check if [".$configHash->{'unirefDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
                  #$wrongLinks++;
                #}
              #}
            #}
          #}
        #}
        ##print "Uniref DB version link...";
        #if (!head($configHash->{'unirefVerLink'}) ) {
          #if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #if (testFTP($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                #my $new_link = check_FTP_diff_sources($configHash->{'unirefVerLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                #if ( $new_link ne ''){
                    #$configHash->{'unirefVerLink'} = $new_link;
                #}else{
                  #print_and_log ("Please check if [".$configHash->{'unirefVerLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
                  #$wrongLinks++;
                #}
              #}
            #}
          #}
        #}
        ##print "Uniprot idmapping link...";
        #if (!head($configHash->{'GODBLink'}) ) {
          #if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
            #if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
              #if (testFTP($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'}) != 1) {
                #my $new_link = check_FTP_diff_sources($configHash->{'GODBLink'}, $configHash->{'uniprotWebUser'}, $configHash->{'uniprotWebPass'});
                #if ( $new_link ne ''){
                    #$configHash->{'GODBLink'} = $new_link;
                #}else{
                  #print_and_log ("Please check if [".$configHash->{'GODBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
                  #$wrongLinks++;
                #}
              #}
            #}
          #}
        #}
        
        ##Links that are not to Uniprot will use double time he Head function and finally the Net::FTP check
        ##NOT UNIPROT
         #if (!head($configHash->{'enzymeDBLink'}) ) {
          #if (!head($configHash->{'enzymeDBLink'})) {
            #if (testFTP($configHash->{'enzymeDBLink'})!= 1) {
              #print_and_log ("Please check if [".$configHash->{'enzymeDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file);#
              #$wrongLinks++;
            #}	
          #}
        #}
        #if (!head($configHash->{'cdDBLink'})) {
          #if (!head($configHash->{'cdDBLink'})) {
            #if (testFTP($configHash->{'cdDBLink'}) != 1) {
              #print_and_log ("Please check if [".$configHash->{'cdDBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
              #$wrongLinks++;
            #}
          #}
        #}
        ##This link always has some more troubles... let's add one more check!
        #if (!head($configHash->{'cdTableLink'})) {
          #if (!head($configHash->{'cdTableLink'})) {
            #if (!head($configHash->{'cdTableLink'})) {
              #if (testFTP($configHash->{'cdTableLink'}) != 1) {
                #print_and_log ("Please check if [".$configHash->{'cdTableLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
                #$wrongLinks++;
              #}
            #}
          #}
        #}
        

        #if (!head($configHash->{'ncRNADBLink'})) {
          #if (!head($configHash->{'ncRNADBLink'})) {
            #if (testFTP($configHash->{'ncRNADBLink'}) != 1) {
              #print_and_log ("Please check if [".$configHash->{'ncRNADBLink'}."] is a correct URL. Annocript will continue...\n",$log_file); 
              #$wrongLinks++;
            #}
          #}
        #}

        #if (!head($configHash->{'pathwaysTableLink'})) {
          #if (!head($configHash->{'pathwaysTableLink'})) {
            #if ( testFTP($configHash->{'pathwaysTableLink'}) != 1 ) {
              #print_and_log ("Please check if [".$configHash->{'pathwaysTableLink'}."] is a correct URL. Annocript will continue...\n",$log_file); #
              #$wrongLinks++;
            #}
          #}
        #}
      #If some check went wrong
      if ($wrongLinks > 0){
        my $text = "Since there is a problem with the connection, Annocript will try different kinds of downloads. Please try on a browser if the links above work".
        " and then press y to continue.\n";
        print_and_log($text,$log_file);#Print on the log file      
       }else{  
				print_and_log ("Everything is ok!\n",$log_file);
		  }
  }
  #This is a simple control that is done even if the database is not created. It is more simple and faster.
  elsif (! (checkLink($configHash->{'ncRNADBLink'}) and checkLink($configHash->{'cdDBLink'}) and checkLink($configHash->{'enzymeDBLink'}) and 
          checkLink($configHash->{'GOTermsLink'}) and checkLink($configHash->{'GODBLink'})
            and checkLink($configHash->{'unirefVerLink'}) and checkLink($configHash->{'uniprotkbVerLink'}) and checkLink($configHash->{'uniprotVerLink'}) 
           # and checkLink($configHash->{'silvaLSULink'}) and checkLink($configHash->{'silvaSSULink'}) ) )
            and checkLink($configHash->{'tremblDBLink'}) and checkLink($configHash->{'swissprotDBLink'}) ) ){
                print_and_log ("Pay attention... you did some error in writing links in the configuration file. Please check it...\n",$log_file);
  }elsif ( ($configHash->{'doBuildOutput'} eq 'YES') or ($configHash->{'doBuildOutput'} eq 'YES' )){
    
	  #The first check is on the permissions. If the user does not have, he cannot create database or perform analysis
	  if (!check_DB_user_permissions($configHash->{'platform'},$configHash->{'host'},
				  $configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'})) {
		 log_and_exit("ERROR: You do not have enough privileges to create the DB. Please ask your adiministrator. Annocript will close...\n",$log_file);
	  }
  }
  
  # print_and_log ("Checking CD db..\n",$log_file);#DEBUGCODE
  #Check if the database chose for the mostly expressed domains is correct
  my $cdName4Expression = $configHash->{'cdName4Expression'};
  #print "cdName4Expression : $cdName4Expression\n";
  if ( !(grep {/$cdName4Expression/ } @cdDBAdmitted) ){
    log_and_exit("Can't find CD database $cdName4Expression. You should "
        ."choose one among:\n ".print_array(\@cdDBAdmitted)."\nAnnocript will exit...\n",$log_file);
  }
      
  #Check on the log file for 
  # print_and_log ("Checking log files..\n",$log_file);#DEBUGCODE
  if ( !($configHash->{'headParsNAValues'} =~ /\.log?/)){
    log_and_exit( $configHash->{'headParsNAValues'}." is not a correct log file name. It must be something like: file.log\n",$log_file);
  }
 
  if ( !($configHash->{'uniprotGenesNotPresent'} =~ /\.log?/)){
    log_and_exit( $configHash->{'uniprotGenesNotPresent'}." is not a correct log file name. It must be something like: file.log\n",$log_file);
  }
  
  #print_and_log ("Checking userid and pwd..\n",$log_file);#DEBUGCODE
  #A check on the user-id and password: it can be of every type of character so 
  #it don't needs a particular control
  if ( !(defined ($configHash->{'mySqlUser'})) or !(defined ($configHash->{'mySqlPass'}) ) ){
    log_and_exit( "ERROR: You need to set an user-id and password for MySQL. Please contact your adiministrator".
    " if you can't. Annocript will exit...\n",$log_file);
  }#else{print_and_log ("MYSQL account: ".$configHash->{'mySqlUser'}." and pwd: ".$configHash->{'mySqlPass'}."\n",$log_file);}#DEBUGCODE
  
  #A check on the parsing file path
  if ( -e $dataFolder.$configHash->{'parsingFilePath'} ){
    log_and_exit( "ERROR : File ".$configHash->{'parsingFilePath'}." already exists.\n",$log_file); #Please change the parsingFilePath in ".$configFilePath."\n";
  }
              
  #print_and_log ("Checking configuration..\n",$log_file);#DEBUGCODE

  #If the user sets the config file to not execute the DB_CREATION, this means that nothing will 
  #be downloaded: the user has already a database.
  my $dbDataFolder; 
    if ( $configHash->{'doDbCreation'} ne "YES"){
        if( !(scalar <$mainDataFolder/*>) ){
           log_and_exit( "ERROR: data folder is empty... You must run Annocript with doDbCreation = YES. Closing",$log_file);
        }else{
          #In the first case Annocript will ask to the user what folder to use as database
          print_and_log ("DB Creation not active: Database will not be downloaded and installed.\n",$log_file);
          
          $configHash->{'dbVersion'} = '';
          
          #If the user wants to use a previous execution of softwares, then the database used can be 
          #only the one used before
          if ( ($configHash->{'doExecutePrograms'} eq "NO") and 
          ( ($configHash->{'doBuildOutput'} eq "YES") or ($configHash->{'extractStatistics'} eq "YES") ) ){
						my $versions_file = $dataFolder."/".$configHash->{'versions_file'};
						print "Checking versions used previously in file $versions_file\n";
						#The database version is present in the version file
						open (VER_F,"<$versions_file")
							or die "Cannot open $versions_file. Maybe programs have not been executed before..\n";
						while (my $row = <VER_F>){
							if ( $row =~ /UniProt version/){
								chomp($row);
								$configHash->{'dbVersion'} = (split(":",$row))[1];
								print_and_log("Version of the database used previously is ".$configHash->{'dbVersion'}.".\n",1,$log_file);
							}
							if ( $row =~ /UniProt database used/){
								chomp($row);
								$configHash->{'dbInUse'} = (split(":",$row))[1];
								print_and_log( "UniProt database used previously is ".$configHash->{'dbInUse'}.".\n",$log_file);
							}
							if ( $row =~ /Domains database used/){
								chomp($row);
								$configHash->{'cdDB'} = (split(":",$row))[1];
								print_and_log( "Domains database used previously is ".$configHash->{'cdDB'}.".\n",$log_file);
							}
							if ( $row =~ /UniRef type used/){
								chomp($row);
								$configHash->{'unirefDB'} = (split(":",$row))[1];
								print_and_log( "UniRef version used is ".$configHash->{'unirefDB'}.".\n",$log_file);
							}
						}
						close(VER_F);
						if ( $configHash->{'dbInUse'} eq "uniref" ){ 
								$configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
						}else{ 
								$configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
						}
						#Check if database used is present
						$configHash->{'database'} = $configHash->{'dbInUse'}.'_'.$configHash->{'dbVersion'};#Gives the name to the database
						log_and_exit("Annocript cannot find the mysql database".$configHash->{'database'}."\n",$log_file) unless ( db_present($configHash->{'database'},$configHash->{'platform'},$configHash->{'host'},
						$configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'}) > 0);
						print_and_log( "MySQL database in use is: ".$configHash->{'database'}."\n",$log_file);
						$dbDataFolder = $mainDataFolder."/".$configHash->{'dbVersion'};
					}elsif ( ($configHash->{'doExecutePrograms'} eq "NO") and 
						($configHash->{'doBuildOutput'} eq "NO") and ($configHash->{'extractStatistics'} eq "NO") ){
							log_and_exit("You are running nothing! Please use the user guide to learn how to use Annocript.\n",$log_file);
					}
          
			 #The user wants to execute the programs, therefore he can choose the database
			if ($configHash->{'dbVersion'} eq ''){
				if ( $db_version ne '' and $db_type ne ''){
					my $dbName = "$db_type\_$db_version";
					#Here we access the mysql system to see if the database is present
					my $dbPresent = db_present($dbName,$configHash->{'platform'},$configHash->{'host'},
					$configHash->{'port'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
					
					if ( $dbPresent > 0){
						#Gives the name to the database
						$configHash->{'database'} = $dbName;	
						$configHash->{'dbInUse'} = $db_type;
						$configHash->{'dbVersion'} = $db_version;
						#The CD database if it is present
						$configHash->{'cdDB'} = "Cdd";#choose_CD_database($dbDataFolder);
						if ( $configHash->{'dbInUse'} eq "uniref" ){ 
								$configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
						}else{ 
								$configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
						}							
						$dbDataFolder = $mainDataFolder."/".$configHash->{'dbVersion'};
						if ($uniref_db_ver eq ''){
							$uniref_db_ver = extract_name($configHash->{'unirefDBLink'},1);
							$configHash->{'unirefDB'} = $uniref_db_ver.".fasta";
						}else{
							$configHash->{'unirefDB'} = $db_type.$uniref_db_ver.".fasta";		
						}
					}							
				}else{
					log_and_exit("WARNING: You must choose an UniProt version of database among those you already have with --db_version. Exiting\n",$log_file);
				}
			}

			#We store the names for all the database. They will not be downloaded but we need them later  
			$configHash->{'tremblDB'} = extract_name($configHash->{'tremblDBLink'},'gz');
			$configHash->{'swissProtDB'} = extract_name($configHash->{'swissprotDBLink'},'gz');
			$configHash->{'cdDBPath'} = $dbDataFolder."/".$configHash->{'cdDB'};			
			$configHash->{'uniprotKBDB'} = "uniprot_kb_".$configHash->{'dbVersion'}.".fasta";#Creates the name for the uniprot kb	
	
				
		 #Here I do a cross check. If the version of the database used does not correspond
		 #with the database chosen to use, it will give an error
		 my $dsn = "dbi:".$configHash->{'platform'}.":".$configHash->{'database'}.":".$configHash->{'host'}.":".$configHash->{'port'};
		 check_db_tables_version($dsn,$configHash->{'mySqlUser'},$configHash->{'mySqlPass'},$dbTablesVersion) 
			or log_and_exit( "The version of the database you selected is not compatible with this version of Annocript ($version) ". 
			" Please choose another or create a completely new database. Exiting...\n",$log_file);  	
				
	     
	     #Then do many checks on the indexed databases and on the choose of a single or many organisms but only
	     #if the user wants to execute analyses, otherwise it does not matter
	     if ($configHash->{'doExecutePrograms'} eq 'YES'){     
				 #Changing the database name if a single organism is to be used, only when using TrEMBL
				 if ( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniprotkb") ){
						#If the file with the sequences exists in the working folder... 
						if (-e $workingFolder.'/'.$configHash->{'blastedOrganism'}){
							my $blOrg = extract_name($configHash->{'blastedOrganism'},0);
							
							$configHash->{'tremblDB'} = "TR_$blOrg.fasta";
							$configHash->{'swissProtDB'} = "SP_$blOrg.fasta";
							#check if there is a file named in same way in the folder. if it is there, the user is alerted
							if (-e $dbDataFolder."/".$configHash->{'tremblDB'}){
								print_and_log(	"The indexed database named ".$configHash->{'tremblDB'}." is present in".
								" $dbDataFolder. Annocript will use it or overwrite if not complete...\n",$log_file);
							}
							if (-e $dbDataFolder."/".$configHash->{'swissProtDB'}){
								print_and_log(	"The indexed database named ".$configHash->{'swissProtDB'}." is present in".
								" $dbDataFolder. Annocript will use it or overwrite if not complete...\n",$log_file);
							}
							copy($workingFolder.'/'.$configHash->{'blastedOrganism'},$dbDataFolder.'/'.$configHash->{'blastedOrganism'});				
						}else{
							log_and_exit("There is no file named ".$configHash->{'blastedOrganism'}." in $workingFolder. Please".
							" check you wrote a correct name in the config file and that the file with organisms names is in your".
							" working folder. If you want to blast against all organisms please change the value with 'all'. \n",$log_file);
						}  
					}elsif( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniref")){
						log_and_exit("When using UniRef all the organisms will be evalutated for the blast. If you want to blast against ".
						" specific organisms you should use TrEMBL.\n",$log_file);
					 }
					
					#Check on the databases. The index file should be present, otherwise the user will be alerted that
					#he should start dbCreation procedure
					if ($configHash->{'dbInUse'} eq "uniprotkb"){
						if ( (indexed_db_present($configHash->{'tremblDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
						 log_and_exit("Database ".$configHash->{'tremblDB'}." is not indexed. Please restart annocript with doDbCreation=YES".
						 "or change the parameter blastedOrganism in your configuration file if you want to use a previously created ".
						 " database of organisms\n",$log_file);
						}				
					}
					
					if ($configHash->{'dbInUse'} eq "uniref"){
						if ( (indexed_db_present($configHash->{'unirefDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
						 #print_and_log("Database ".$configHash->{'uniref'}." is not indexed. Annocript will do it!\n",$log_file);
						 log_and_exit("Database ".$configHash->{'unirefDB'}." is not indexed. Please restart annocript with doDbCreation=YES\n",$log_file);
						}				
					}
					
					if ((indexed_db_present($configHash->{'swissProtDB'}, "blastIndexedExts",$dbDataFolder) == 0) ){
					 log_and_exit("Database ".$configHash->{'swissProtDB'}." is not indexed. Please restart annocript with doDbCreation=YES\n",$log_file);
					}
	   }#END IF the user want to execute programs
	 }#END If the data folder is empty
  }#END IF the user doesn't want to create the database
  else{
    #The database has to be create. The most recent version will be downloaded
    #We need a value for db_type because db_version will be the last one
		if ( $db_type eq ''){
		 log_and_exit("Please use the parameter --db_type to select which database you want to install [uniref or uniprotkb (for trembl)].".
		 " Consider that if you want to blast against specific species you need uniprotkb.\n",$log_file);				
		}
		$configHash->{'dbInUse'} = $db_type;
		if ( $configHash->{'dbInUse'} eq "uniref" ){ 
			$configHash->{'versionDBLink'} = $configHash->{'unirefVerLink'};
      $configHash->{'unirefDB'} = extract_name($configHash->{'unirefDBLink'},'gz');				
		}else{ 
			$configHash->{'versionDBLink'} = $configHash->{'uniprotVerLink'};
      $configHash->{'tremblDB'} = extract_name($configHash->{'tremblDBLink'},'gz');				
		}	
		
		#First case: the user wants to install a new database
		if ( $db_version eq ''){
		 print_and_log("You did not choose the version of the db using --db_version. The latest will be used...\n",$log_file);	
		 $db_version = getDBVersion($configHash->{'versionDBLink'},$mainDataFolder,$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'}, $configHash->{'dbInUse'});
			print_and_log ("\nCurrent UniProt db version will be installed: ".$db_version."\n",$log_file);
			$configHash->{'dbVersion'} = $db_version;		
		}
		#else he wants to update an existing one
		else{
			if ( -d $mainDataFolder."/".$db_version){
				print_and_log ("\nDatabase $db_version will be overwritten\n",$log_file);
				$configHash->{'dbVersion'} = $db_version;						
			}else{
				log_and_exit("ERROR: the database $db_version do not exist. Choose an existing one from $mainDataFolder. Exiting\n",$log_file);
			}
		}
		my $dbName = "$db_type\_$db_version";	

    #The name of the SwissProt db
    $configHash->{'swissProtDB'} = extract_name($configHash->{'swissprotDBLink'},'gz');		
	
    #..and the folder.
    $dbDataFolder = $mainDataFolder."/".$db_version;#Create the name for the folder..
    $configHash->{'database'} = $configHash->{'dbInUse'}.'_'.$db_version;#Give the name to the database
    print_and_log ("Folder for DB will be: ".$dbDataFolder." and database name: ".$configHash->{'database'}."\n",$log_file);
    $tempTableFolder =  $dbDataFolder."/".$tempTableFolder;#Creating the tables folder
    $configHash->{'tempTableFolder'}  = $tempTableFolder;#Creating a global variable to share among modules
    
    #Creating the db folder
    unless(-d $dbDataFolder){
      print_and_log ($dbDataFolder." doesn't exists. Creating folder $dbDataFolder...\n",$log_file);
      mkdir $dbDataFolder or log_and_exit( "ERROR : Can't create $dbDataFolder. Check permissions!\n",$log_file);
    }else{#If the folder is there, check the UniRef db he wants to install
			if ($configHash->{'dbInUse'} eq "uniref"){
				#Check if there is a uniref database in the db folder, if something is there the version of uniref is compared
				#with the uniref link and if it is different the user is alerted
				my @ufPresent = check_uniref_presence($dbDataFolder);
				if (scalar(@ufPresent) > 0){
					my $fromLink = extract_name($configHash->{'unirefDB'},1);
					if( !(grep  {/$fromLink/} $ufPresent[0]) ){
					log_and_exit( "The type of UniRef you want to download is not the same as the one you downloaded in a previous execution".
					" of Annocript. Annocript can manage only one type of UniRef (50,90 or 100) for each versions of the Uniprot db.".
					"Check the unirefDBLink and unirefVerLink in CONFIGURATION/config_annocript.txt file.\n",$log_file);
					}
				}
			}
			print_and_log ("A folder $dbDataFolder already exists. Annocript will continue...\n ",$log_file);
		}

			
 #Creates a folder where all the tables will be put
	unless(-d $tempTableFolder){
		print_and_log ($tempTableFolder." doesn't exists. Creating folder $tempTableFolder...\n",$log_file);
		mkdir $tempTableFolder or log_and_exit( "ERROR : Can't create $tempTableFolder. Check permissions!\n",$log_file);
	}else{print_and_log ("A folder $tempTableFolder already exists. Annocript will continue...\n ",$log_file);}
	
	$configHash->{'uniprotKBDB'} = "uniprot_kb_".$db_version.".fasta";#Create the name for the uniprot kb
	
	#Move the file with the database version to the folder just created
	move($mainDataFolder."/".(extract_name($configHash->{'versionDBLink'},0)),$dbDataFolder); 
		
	#Give the  name to the CDD database extracting from the link
	$configHash->{'cdDB'} = extract_name($configHash->{'cdDBLink'},"1");
	#Sometimes in uniprot they put the name of the zip file with LE (Little Endian) and BE (Big Endian)
	$configHash->{'cdDB'} = clean_cd_name($configHash->{'cdDB'});
	 
	$configHash->{'cdDBPath'} = "$dbDataFolder/".$configHash->{'cdDB'};
	
	print_and_log ("Conserved domains database: ".$configHash->{'cdDB'}."\n",$log_file);
	
	if ($configHash->{'eraseFiles'} eq "YES"){
			print_and_log ("Since you selected eraseFiles=YES into the configuration file,  ".
						" Annocript will remov the files downloaded for the database generation and BLAST execution.\n",$log_file);
	}
	
	#Here Annocript decides if to start the DB creation using more memory and doing fast or not
	#It depends on how much memory the system has but it can be choose inside the annocript configuration file
	my $ram = ram_memory();
	if ( $ram/1024/1024 > $configHash->{'bigMemory'} ){
		$configHash->{'useHugeRAM'} = 'YES';
		print_and_log("You have enough memory to create more rapidly the db. Annocript will use a lot of memory now! :-D\n",$log_file);
	}else{
		$configHash->{'useHugeRAM'} = 'NO';
		print_and_log("You don't have enough memory to be more rapid. Annocript will be little slow but won't use lot of memory! :-)\n",$log_file);
	}
      
     #Changing the database name if selected organims will be used
	 if ( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniprotkb") ){
	  #Search the file with organisms	
	  if (-e $workingFolder.'/'.$configHash->{'blastedOrganism'}){
		my $blOrg = extract_name($configHash->{'blastedOrganism'},0);
		$configHash->{'tremblDB'} = "TR_$blOrg.fasta";
		$configHash->{'swissProtDB'} = "SP_$blOrg.fasta";
		#check if there is a file named in same way in the folder. if it is there, the user is alerted		
	    if (-e $dbDataFolder."/".$configHash->{'tremblDB'}){
			print_and_log("An existing database ".$configHash->{'tremblDB'}." is present in the folder ".
			"$dbDataFolder. Please change the name of the file of organisms or remove this file.\n" ,$log_file );
			exit;
		}
		if (-e $dbDataFolder."/".$configHash->{'swissProtDB'}){
			print_and_log("An existing database ".$configHash->{'swissProtDB'}." is present in the folder ".
			"$dbDataFolder. Please change the name of the file of organisms or remove this file.\n"  ,$log_file);
			exit;
		}		
		#Shift the file with the organisms to blast inside the db folder
		copy($workingFolder.'/'.$configHash->{'blastedOrganism'},$dbDataFolder.'/'.$configHash->{'blastedOrganism'});
				
	  }else {
		log_and_exit("There is no file named ".$configHash->{'blastedOrganism'}." in $workingFolder. Please".
		 " check you wrote a correct name in the config file and that the file with organisms names is in your".
		 " working folder. \n",$log_file);
	  } 
	 }elsif( ($configHash->{'blastedOrganism'} ne 'all') and ($configHash->{'dbInUse'} eq "uniref")){
			log_and_exit("When using UniRef all the organisms will be evaluated for the blast. If you want to blast against ".
				" specific organisms you should use Trembl. Exiting...\n",$log_file);
	  }
			    
  }#END OF CODE WHEN THE DATABASE IS TO BE GENERATED
  	
  if (defined $dbDataFolder){
    $configHash->{'dbDataFolder'} = $dbDataFolder;
  }
  
  #A check on the database name...
  uniprot_db_selection($configHash,$log_file);
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
  
  #annoPrint (Dumper\$configHash;
	#CONTROL THE PATH TO DATABASES  
	if ( $configHash->{'ncDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name".$configHash->{'ncDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n",$log_file);
	}   
	if (  $configHash->{'swissProtDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name".$configHash->{'swissProtDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n",$log_file);
	}
	if ( $configHash->{'tremblDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name".$configHash->{'tremblDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-\.]. Program will exit...\n",$log_file);
	}
             
             
  #Check if the links wrote for the programs are valid
  check_programs_presence($configHash,\@blastProgramsInUse,$log_file);     
  
  #Check the blast input parameters using the final chars near the parameter
  #X is used for BLASTX and BLASTP, N for BLASTN and TBLASTN, RPS for RPSTBLASTN and RPSBLAST
  my @blastIds = ('X','N','RPS'); 
          
  foreach my $blastId (@blastIds){
    my $arg = 'word_size'.$blastId;
			if (defined($configHash->{$arg})){  
				 if (!correct_type($configHash->{$arg},"positiveint")){
          log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
         }  
      }
      $arg = 'evalue'.$blastId;
      if (defined($configHash->{$arg})){  
        if (!correct_type($configHash->{$arg},"real")){
          log_and_exit( $configHash->{$arg}." is not a real. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}  
			}
			$arg = 'num_descriptions'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}  
			}
			$arg = 'num_alignments'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{ $arg},"positiveint")){
					log_and_exit( $configHash->{ $arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}  
			}
			$arg = 'max_target_seqs'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}  
			}
      $arg = 'num_threads'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}else{
          my $num_processors = num_processors();
       
          if ($configHash->{$arg} > $num_processors ) {
              log_and_exit ("ERROR: Please change the paramater $arg in config file. Keep in mind that your processors are: $num_processors.\n",$log_file);
          }
        }  
			}
      $arg = 'threshold'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}  
			}
      $arg = 'outfmt'.$blastId;
			if (defined($configHash->{$arg})){  
				if (!correct_type($configHash->{$arg},"positiveint")){
					log_and_exit( $configHash->{$arg}." is not a positive integer. Please change variable $arg and re-execute Annocript. Exiting...\n",$log_file);
				}else{
            
        }  
			}
		}
    
    #Number of processes to use for the parallel execution 
    #This checks if your choice is good depending on the number of processes used      
		if (defined($configHash->{'threads4Parallel'})){  
			if (!correct_type($configHash->{'threads4Parallel'},"positiveint")){
				log_and_exit( $configHash->{'threads4Parallel'}." is not a positive integer. Please change threads4Parallel and re-execute Annocript. Exiting...\n",$log_file);
			}else { 
        my $num_processors = num_processors();
       
        if ($configHash->{'threads4Parallel'} > $num_processors ) {
            log_and_exit ("ERROR: Please change the paramater threads4Parallel in config file. Keep in mind that your processors are: $num_processors.\n",$log_file);
        }      
      }  
		}
    	 
		# Chekcing correctness of lncRNAPrediction's parameters
		if (defined($configHash->{'NCSeqLength'})){  
			if (!correct_type($configHash->{'NCSeqLength'},"positiveint")){
				log_and_exit( $configHash->{'NCSeqLength'}." is not a positive integer. Please change NCSeqLength and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		if (defined($configHash->{'NCORFLength'})){  
			if (!correct_type($configHash->{'NCORFLength'},"positiveint")){
				log_and_exit( $configHash->{'NCORFLength'}." is not a positive integer. Please change NCORFLength and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		if (defined($configHash->{'NCThresh'})){  
			if (!correct_type($configHash->{'NCThresh'},"real")){
				 log_and_exit( $configHash->{'NCThresh'}." is not a real. Please change NCThresh and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		
    if (defined($configHash->{'max_port_seq_length'})){  
			if (!correct_type($configHash->{'max_port_seq_length'},"positiveint")){
				log_and_exit( $configHash->{'max_port_seq_length'}." is not a positive integer. Please change max_port_seq_length and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
    if (defined($configHash->{'min_port_seq_length'})){  
			if (!correct_type($configHash->{'min_port_seq_length'},"positiveint")){
				log_and_exit( $configHash->{'min_port_seq_length'}." is not a positive integer. Please change min_port_seq_length and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		
		#Checking the number of tentatives to execute commands
	  if (defined($configHash->{'max_tentatives'})){ 
			if (!correct_type($configHash->{'max_tentatives'},"positiveint")){
				log_and_exit( $configHash->{'max_tentatives'}." is not a positive integer. Please change max_tentatives and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		
    #Checking parameters to use for dna2pep
    if (defined($configHash->{'d2pMode'})){  
			if ( ($configHash->{'d2pMode'} ne "strict") and ($configHash->{'d2pMode'} ne "any") and ($configHash->{'d2pMode'} ne "none") ){
				 log_and_exit(" d2pMode should be one among: strict, any, none. And it is: ". $configHash->{'d2pMode'}.". Please change it and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}    
      
    #Checking the maximum allowed length for the rows in the fasta file
    if (defined($configHash->{'maxFastaSeqLen'})){  
			if (!correct_type($configHash->{'maxFastaSeqLen'},"positiveint")){
				log_and_exit( $configHash->{'maxFastaSeqLen'}." is not a positive integer. Please change maxFastaSeqLen and re-execute Annocript. Exiting...\n",$log_file);
			}  
		}
		
    #Annocript Checks if there is a file of counts and if it is good
    #check_counts_file();
    #IF ONE DAY WE'LL DECIDE TO USE THE COUNTS WE CAN START FROM HERE... TILL THAT...
    $configHash->{'countOutPresent'} = 'NO';
}



=head2 checkPROGRAMS_OUT_PARSER

 Title   : checkPROGRAMS_OUT_PARSER
 Usage   : checkPROGRAMS_OUT_PARSER(  
								- local folder = the folder where actually is the checkPROGRAMS_OUT_PARSER script
                               );

 Function: Checks the config hash parameters regarding the ouput creation: how they are written, if they are present and so on..
 Returns : nothing

=cut
sub checkPROGRAMS_OUT_PARSER{
  
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
  $configHash->{'destroy_gffdb'} = detect_gff_db($configHash,$log_file);
		
  #A check of the file extensions...
  if ( !($configHash->{'blastxGFF3FileName'} =~ /\.gff?/)){
    log_and_exit($configHash->{'blastxGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n",$log_file);
  }  
  if ( !($configHash->{'blastxTRGFF3FileName'} =~ /\.gff?/)){
    log_and_exit($configHash->{'blastxTRGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n",$log_file);
  } 
  if ( !($configHash->{'blastxSPGFF3FileName'} =~ /\.gff?/)){
    log_and_exit( $configHash->{'blastxSPGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n",$log_file);
  } 
  if ( !($configHash->{'rpstblastnGFF3FileName'} =~ /\.gff?/)){
    log_and_exit( $configHash->{'rpstblastnGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n",$log_file);
  } 
  if ( !($configHash->{'blastnGFF3FileName'} =~ /\.gff?/)){
    log_and_exit( $configHash->{'blastnGFF3FileName'}." is not a correct GFF file name. It must be something like: file.gff\n",$log_file);
  }
    
  if ( !($configHash->{'ORFFastaFileName'} =~ /(\.fasta?)|(\.fa?)/)){
    log_and_exit( $configHash->{'ORFFastaFileName'}." is not a correct FASTA file name. It must be something like: file.fasta\n",$log_file);
  } 
  if ( !($configHash->{'fastaSeqs'} =~ /(\.fasta?)|(\.fa?)/)){
    log_and_exit( $configHash->{'fastaSeqs'}." is not a correct FASTA file name. It must be something like: file.fasta\n",$log_file);
  }
  if ( !($configHash->{'outFileName'} =~ /\.txt?/)){
    log_and_exit( $configHash->{'outFileName'}." is not a correct txt file name. We prefer to have a txt file for all results: output.txt\n",$log_file);
  }
  #A control on the type of the output of blast programs
  if ( !($configHash->{'blastOutType'} =~ /blast/)){
    log_and_exit( $configHash->{'blastOutType'}." the blast out type supported is blast, keep in mind and modify the configuration file!\n",$log_file);
  }
  
    #Checking the variable plotType. It can have only one of two values
  if ( !($configHash->{'plotType'} =~ /(piechart?)|(barplot?)/)){
    log_and_exit( "You can choose only one between [piechart, barplot] as a value for plotType in the configuration file. ".$configHash->{'plotType'}." is not allowed!\n",$log_file);
  } 
  
  #Checking the variable goTermsAss. It can have only one of three values
  if ( !($configHash->{'goTermsAss'} =~ /(proteins?)|(domains?)|(both?)/)){
    log_and_exit( "You can choose only one between [proteins, domains] as a value for goTermAss in the configuration file. ".$configHash->{'goTermsAss'}." is not allowed!\n",$log_file);
  } 
}
####################################


=head2 set_global_var_and_confighash

 Title   : set_global_var_and_confighash
 Usage   : set_global_var_and_confighash(   );

 Function:  Just sets important variables and imports the parameters from the configuration 	
					file into the hash
 
 Returns : nothing

=cut
sub set_global_var_and_confighash{
	my $config_user = shift;
	
	print_and_log($program_title,$log_file);
	##Set the folders name using inizialize_folders
	print_and_log($configClosure,$log_file);	
	($workingFolder,$program_folder) = initialize_folders($foldersFile);
	$config_folder = $program_folder."/".$config_folder;
	$program_config = $config_folder."/".$program_config;
	$variables_file = $config_folder."/".$variables_file;
	

	##Set important variables
	$mainDataFolder = $workingFolder."/".$mainDataFolder;
	$analysesFolder = $workingFolder."/".$analysesFolder;
   #print_and_log("analysesFolder $analysesFolder  working_folder  $workingFolder program_folder $program_folder \n",$log_file);#DEBUGCODE
	
	#Check if config files are ok (1 and 2 indicates the line into the variables file
	print_and_log("Checking $program_config \n",$log_file);	
	#print_and_log("using line 1 of $variables_file...\n",$log_file);	#DEBUGCODE
	checkConfigVariables($program_config,$variables_file,1);
	print_and_log("Checking $config_user \n",$log_file);	 
	#print_and_log("using line 2 of $variables_file...\n",$log_file);	#DEBUGCODE
	checkConfigVariables($config_user,$variables_file,2);
	
	#loads the parameters in the configuration hash as from the configuration files
	print_and_log("Loading the parameters for the run as from $program_config and $config_user...\n",$log_file);	
	configFile2Hash($program_config,\$configHash);
	configFile2Hash($config_user,\$configHash);
	print_and_log("DONE!\n",$log_file);	#DEBUGCODE
	
  $configHash->{'dataFolder'} = $dataFolder;
  $configHash->{'gffFolder'} = $gffFolder;
  $configHash->{'outFolder'} = $outFolder;
  $configHash->{'statsFolder'} = $statsFolder;
  $configHash->{'logFolder'} = $logFolder;
	#Storing the session folder in the hash  
	$configHash->{'analysesFolder'} = $analysesFolder;

                  	  
	#Same some variable in the config hash
	$configHash->{'work_fold'} = $workingFolder;
	$configHash->{'lib_folder'} = $program_folder."/".$lib_folder;
	$configHash->{'usefulFolder'} = $program_folder."/".$configHash->{'usefulFolder'};
	$configHash->{'guideFolder'} = $program_folder."/".$configHash->{'guideFolder'};
	$configHash->{'configFolder'} = $program_folder."/".$configHash->{'configFolder'};
	
	#print_and_log("Making folders data and analyses!\n",$log_file);	#DEBUGCODE
  #Creates the jobs folder if it does not exists
  unless(-d $analysesFolder){
    print_and_log ($analysesFolder." doesn't exists. Creating folder...\n",$log_file);
    mkdir $analysesFolder or log_and_exit("ERROR: can't create folder $analysesFolder\n ",$log_file);
  }
  #Creates the main 'data' folder if it does not exists  
  unless(-d $mainDataFolder){
    print_and_log ($mainDataFolder." doesn't exists. Creating folder...\n",$log_file);
    mkdir $mainDataFolder or log_and_exit( "ERROR: can't create folder $mainDataFolder\n ",$log_file);
  }
	
	#print Dumper\$configHash;#DEBUGCODE
}

=head2 check_input_parameters

 Title   : check_input_parameters
 Usage   : check_input_parameters(  );

 Function:  Identifies the analysis that the user is going to use
 
 Returns : 1 if the user choose a new session else 0

=cut	
sub check_input_parameters{
	
	#If the user did not give any analysis name should mean that he only wants
	#to build the database or make some utils
	if ($analysis){
		$analysisFolder = $analysesFolder."/".$analysis;
		$configHash->{'sessionFolder'} = $analysisFolder;
		#The analysis already exists  
		if (-d $analysisFolder){
			if (valid_analysis($configHash,$analysisFolder) < 0){
			 log_and_exit( "The $analysis is not valid. Please remove it.\n",$log_file);
			}
			 print_and_log (" Using analysis $analysis...\n",$log_file); 
		}else{
			
			#Check that the user does not insert strange chars in the name
			if ( $analysis !~ /^[A-za-z0-9\_\-]*$/i ){
				log_and_exit( "The name $analysis is not permitted. \nChoose another using [A-za-z0-9\_\-]:",$log_file);
			}
			
			#Creates the folder
			mkdir $analysisFolder 
				or log_and_exit("ERROR: Unable to create the directory $analysisFolder please check your permissions.Exiting..\n",$log_file);
			print_and_log ("A new analysis ($analysis) has been initialized with folder: $analysisFolder!\n",$log_file); 
			$newSession = 1;
		}		
	}else{
			log_and_exit( "Please use the parameter --analysis to select an analysis with which you start. Exiting..\n",$log_file);
	}
	
	#The paramater $db_type must be always given
	if ( $db_type eq ''){
		log_and_exit("The parameter --db_type is mandatory either when you create a database or when you run an analysis.".
		" Please choose one between uniref and uniprotkb (trembl). Exiting..\n",$log_file);
	}elsif ( ! ($db_type eq 'uniref' or $db_type eq 'uniprotkb' ) ){
		log_and_exit("Allowed values for --db_type are uniref and uniprotkb (trembl). Exiting..\n",$log_file);		
	}	
	#Checking db_version
	if ( $db_version eq ''){
		print_and_log ("A new database will be created by downlading current version from UniProt...\n",$log_file);
	}else{
		print_and_log ("The database $db_version will be reinstalled\n",$log_file); 	
	}
	#If no analysis name is used then is fine to just create the database
	#but if the user wants also to execute analyses then he needs to put the analysisname
	#else{
		#if ( $configHash->{'doDbCreation'} eq 'YES'){
			#if ( ($configHash->{'doExecutePrograms'} eq "YES") or 
          #($configHash->{'doBuildOutput'} eq "YES") or ($configHash->{'extractStatistics'} eq "YES") ){
				#log_and_exit("You cannot run any analysis without the analysis name (use --analysis). Or if you want just".
				#" to create a new database set only doDbCreation=YES. Exiting..\n",$log_file);		
			#}else{
				#if ( $db_type eq ''){
					#log_and_exit("If you want to create a database you should specify if uniref or trembl with --db_type command. Exiting..\n",$log_file);
				#}
				##else{
					##if ( $db_version eq ''){
						##print_and_log ("A new database will be created by downlading current version from UniProt...\n",$log_file);
					##}else{
						##print_and_log ("The database $db_version will be reinstalled\n",$log_file); 	
					##}
			  ##}
			#}
		#}else{
				#print_and_log ("You need to use an analysis name to start a new analysis (use --analysis). Or if you want just".
				#" to create  a new database set doDbCreation=YES\n",$log_file); 
		#}
	#}

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



	my $howToUse = "Use with: perl $program_name \n\n".
	"You may want to create a database and run an analysis or you may want to just run an analysis or just create a database".
	"When you insert an analysis name (with --analysis) if it is already present an old analysis will be re-executed or completed".
	"-c|--user_config: The configuration file with all the parameters needed for the execution. It should have been copied ".
									"into the working folder during the installation.\n".
	"-a|--analysis: The analysis to consider. You can choose one among those in the working directory or start a new one.".
	"-db|--db_version: The version of the database to use. Must be one among those into the DATA folder. When creating a new db".
					" the current version will be installed.\n".					
	"-db_t|--db_type: the db type in UniProt i.e. uniref or uniprotkb". 
	"-uf_db_v|--uniref_db_ver: The version of the UniRef database to use. Must be one among those into the DATA folder. When creating a new db".
					" the default version into the program configuration file will be used.\n".
	"-fr|--force_run: use this parameter whenever you want to force an execution or overwrite an analysis\n".
	"-u|--utils: enters into a menù where you can do some utilities\n".
	"-lic|--licence: shows GNU GPL Open Source warranty and conditions.\n\n";

  #  Parse options
  GetOptions(
           "help" =>        \$HELP,
           "version" =>      \$VERSION,
           "license" =>      \$license,
           "c|user_config=s" => \$user_config,   #It's mandatory
           "a|analysis=s" => \$analysis,   #It's mandatory
           "db|db_version=s" => \$db_version,   #It's mandatory
           "uf_db_v|uniref_db_ver=s" => \$uniref_db_ver, 
           "db_t|db_type=s" => \$db_type,   #It's mandatory
           "u|utils" => \$utils,   #It's mandatory
           
           "fr|force_run" => \$force_run
            );

  #Print a little help
  if ( $HELP ){
    print $howToUse;
    #pod2usage(1);
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "Version: $program_name $version \n";
    exit;
  }  
  
  ## If no arguments were given, then allow STDIN to be used only
  ## if it's not connected to a terminal (otherwise print_and_log (usage)
 # pod2usage("$0: Config file is missing...")  if ( $user_config eq '');
   
  #Exit if no config file is inputed  
  if (! (-e $user_config) ) {
		print $howToUse;
		print "Config file $user_config does not exists. Please write a correct path.\n";
    exit ;
  }

}
