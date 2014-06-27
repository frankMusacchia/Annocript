#######################db_creator.pm###########################################################
# db_creator - A module that creates a database with all informations about UNIPROT sequences. #
################################################################################################
#BE CAREFUL!!!! DON'T FORGET TO UPDATE THE VARIABLE (dbTablesVersion) EACH TIME THE DB STRUCTURE WILL CHANGE!!!

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


package DB_CREATION::db_creator;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( execute_db_creator getKBVersion download_file uncompressFile dl_and_extract extract_name delete_file );

use strict;
use warnings;

use Data::Dumper;#To print our hashes
use File::Copy;#To copy files
use Getopt::Long;
#use Pod::Usage;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Error qw(:try);#To check for error
use Archive::Extract;#To extract every type of archive
use Time::HiRes qw( time );#To keep time
use Switch;#Use of the Switch construct
use List::MoreUtils;

use DBI;#DB library   
use LWP::UserAgent;#Ftp connection library

#Annocript libraries
use PROGRAMS_EXEC::execute_programs qw( check_presence run_blast_indexing check_makeblastdb_out);
use USEFUL::utilities qw(dl_and_extract download_file extract_name delete_file 
                  print_array extract_file my_extract_file try_exec_command deallocate_hash
                  indexing_files_present delete_directory getUnirefType join_files_with_cat);


my $info; #An Hash with a collection of informations..
my $configHash;
 
#Parameters for database
my $platform;
my $host; ## MySQL Host
my $port;
my $dbargs; 
my $dsn;
my $database;
my $mySqlUser;
my $mySqlPass;

#Folders
my $sessionFolder;
my $dataFolder;
my $logFolder;
my $outFolder;
my $dbDataFolder;
my $makeblastdbPath ; 
my $tempTableFolder;#Folder to keep the tables while they are being generated

#Variables to identify the database used. Till now they are for UNIPROT_kb and UNIREFIf you want to add another please
#create SQL code and add here
my $dbDownloaded;
my $unirefSqlScript = 'uniref_db.sql';
my $uniprotSqlScript = 'uniprot_db.sql';


#We need these variables to check if the files of the indexed database are present 
#These are for Swiss-Prot, Trembl, Uniref                           
my @blastIndexedExts = qw(.psq .pin .phr);  

#These are for ncRNA.txt
my @ncRNAExts = qw(.nsq .nin .nhr);  
    
    
    
############################################################################### CODE    

=head2 execute_db_creator

 Title   : execute_db_creator
 Usage   : execute_db_creator( -configHash => the hash used for store the variables used by the module
			       );

 Function: main function of the module. It execute all the subroutines of the modules in the correct order
          to build the database and download and indexing the UNIPROT and UNIREF databases .

 Returns : nothing

=cut
sub execute_db_creator{
  
	#This is an hash table with all the elements read in the config.txt file
	$configHash = shift;
	$sessionFolder = $configHash->{'sessionFolder'};#The session of the user
	$dataFolder = $sessionFolder."/".$configHash->{'dataFolder'};#A folder for data of user
	$outFolder = $sessionFolder."/".$configHash->{'outFolder'};#A folder with outputs from each program
	$dbDataFolder = $configHash->{'dbDataFolder'};#The folder with all the databases formed
  $makeblastdbPath = $configHash->{'makeblastdbPath'};
  $tempTableFolder = $configHash->{'tempTableFolder'};
  
  
  my $dbCreationScript; 
  $dbDownloaded = $configHash->{'dbInUse'};
  
  #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
	#Variables for database
	$database = $configHash->{'database'};
	$mySqlUser = $configHash->{'mySqlUser'};
  $mySqlPass = $configHash->{'mySqlPass'};
	$platform = $configHash->{'platform'};
  $host = $configHash->{'host'};## MySQL Host
  $port = $configHash->{'port'};
           
  $dbargs = {AutoCommit => 0,
                  PrintError => 1};
  $dsn = "dbi:$platform:$database:$host:$port";

	#print "Info written on the config file given in input: \n";
	#print Dumper $configHash;#DEBUGCODE

  print "\n(".scalar(localtime).") >>>>>>>>>>>>>DOWNLOAD OF FILES NECESSARY FOR CD SEARCH AND RIBOSOMIAL<<<<<<<<<<<<<<\n";#LOG
  getRfamAndCdd();
  $partDuration = time-$partTime;
  print "Finished Downloading CD and RFAM DBs in ".$partDuration." seconds\n\n";
  $partTime = time;
		
  if ($dbDownloaded eq 'uniprotkb'){
    $dbCreationScript = $uniprotSqlScript;  
  }
  
  if ($dbDownloaded eq 'uniref'){
    $dbCreationScript = $unirefSqlScript;    
  } 


  #Create DB and tables with SQL query with PERL
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>CREATION OF DATABASE AND ITS TABLES<<<<<<<<<<<<<<\n";#LOG
  if ($configHash->{'dbExistent'} eq 'NO'){
    createDBAndTables($dbCreationScript);
    #print "\n....DONE!\n";
    $partDuration = time-$partTime;
    print "Finished Creating DB and tables in ".$partDuration." seconds\n\n";
    $partTime = time;
  }else{ print "Annocript will not create tables since they exist.\n";}
  
  if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){

   if (is_table_in_db('uniprotkb') == -1){ 
      print "\n(".scalar(localtime).") >>>>>>>>>>>>>CREATION OF THE UNIPROT KNOWLEDGBASE TABLES<<<<<<<<<<<<<<\n";#LOG
      createUniprotKBDB($configHash->{'swissprotDBLink'}, $configHash->{'tremblDBLink'}, $configHash->{'versionDBLink'},$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'},$configHash->{'parsingFilePath'}); #serve per creare UNIPROT
      print "\n....DONE!\n";
      $partDuration = time-$partTime;
      print "Finished Creating Uniprot_kb table in ".$partDuration." seconds\n\n";
      $partTime = time;
    }else { print "uniprotkb tables already present. Going to the next..\n";}
  }
  

  if ( $configHash->{'dbInUse'} eq 'uniref' ){
    print "\n(".scalar(localtime).") >>>>>>>>>>>>>CREATION OF THE UNIREF AND SWISS-PROT TABLES<<<<<<<<<<<<<<\n";#LOG
    print "Annocript will insert both Swiss-Prot informations and Uniref ones.\n";
   
    if (is_table_in_db('uniprotkb') == -1){ 
      print "Creating Swiss-Prot tables...\n";
      createUniprotKBDB($configHash->{'swissprotDBLink'}, undef, $configHash->{'versionDBLink'},$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'},$configHash->{'parsingFilePath'}); #serve per creare UNIPROT
      print "\n....DONE!\n";
      $partDuration = time-$partTime;
      print "Finished Creating Swissprot table in ".$partDuration." seconds\n\n";
      $partTime = time;
    }else { print "Swiss-Prot tables already present. Going to the next..\n";}
    

   if (is_table_in_db('uniref') == -1){ 
      print "Creating Uniref tables...\n";
      createUnirefDB($configHash->{'unirefDBLink'}, $configHash->{'versionDBLink'},$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'},$configHash->{'parsingFilePath'}); #serve per creare UNIPROT
      print "\n....DONE!\n";
      
      $partDuration = time-$partTime;
      print "Finished Creating Uniref table in ".$partDuration." seconds\n\n";
      $partTime = time;
    }else { print "Uniref tables already present. Going to the next..\n";}
  }
  
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>MAPPING OF PROTEINS WITH GO-TERMS, ENZYME TERMS AND PATHWAYS<<<<<<<<<<<<<<\n";#LOG
  if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){  
    
   
   #Check if the table is present and it all go, pathway and enzyme are present or everything should be created again
    if ( is_table_in_db('uniprotkbMappingsTable') == -1 ){
      my $go_enz_create = 1;#A flag to decide if only the hash has to be created or also the table in db
      createMappingTable('uniprotkb','uniprotkbMappingsTable',$go_enz_create);#Inserts all elements related to a mapping of a protein to an identifier
     
      $partDuration = time-$partTime;
      print "Finished the mapping in ".$partDuration." seconds\n\n";
      $partTime = time;
    }else { print "uniprotkbMappingsTable table already present. Going to the next..\n";}
  }
  
  if ( $configHash->{'dbInUse'} eq 'uniref' ){ 
    print "A double mapping will start to fill two mappings tables (Swiss-Prot and Uniref)...\n";
    
    my $go_create = 0;#Initialization. if tables are already  present it will not create.
    #Check if the table is present and it all go, pathway and enzyme are present or everything should be created again
    if ( is_table_in_db('uniprotkbMappingsTable') == -1 ){
      $go_create = 1; #Now the go table has to be create. Later for uniref is useless to create it again
      createMappingTable('uniprotkb','uniprotkbMappingsTable',$go_create);#Inserts all elements related to a mapping of a protein to an identifier
      $go_create = 0;
    }else { print "uniprotkbMappingsTable table already present. Going to the next..\n";}   
        
    #Check if unirefMappingTable is present. 
    if ( is_table_in_db('unirefMappingsTable') == -1 ){
      createMappingTable('uniref','unirefMappingsTable',$go_create);#Inserts all elements related to a mapping of a protein to an identifier
    }else { print "unirefMappingsTable table already present. Going to the next..\n";}   
    
    $partDuration = time-$partTime;
    print "Finished the mapping in ".$partDuration." seconds\n\n";
    $partTime = time;
  }
  
  ##Creating the CD table
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>CREATING A TABLE FOR CONSERVED DOMAINS INFORMATIONS<<<<<<<<<<<<<<\n";#LOG
  #Check if the table is present
  my $cdTable = "cdTable";  
  if (is_table_in_db($cdTable) == -1){
    createCDTable($cdTable,$configHash->{'cdTableLink'},$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'});   
  }else { print "$cdTable table already present. Going to the next..\n";}   
  
  $partDuration = time-$partTime;
  print "Finished the Conserved Domains table in ".$partDuration." seconds\n\n";
  $partTime = time;
 
    
  #If the user is working with an existing database then we need a check of the indexed databases. If everything is ok
  #Annocript can go, else they should be downloaded and installed.
  
  
  #Indexing for BLAST programs
  #Check programs existence
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>INDEXING THE FASTA DATABASES<<<<<<<<<<<<<<\n";#LOG
  print "The following databases are required to execute BLAST programs...\n";
  if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){
    print "Swiss-Prot indexed database:\n";
    if ( indexed_db_present(extract_name($configHash->{'swissProtDB'},0),\@blastIndexedExts) == 0){
      print "Indexing ".$configHash->{'swissProtDB'}."\n";
      #If it is defined the link to swiss-prot then download and extract it
      if (defined $configHash->{'swissprotDBLink'} ){
        my $swissprotDBNameCom = extract_name($configHash->{'swissprotDBLink'}, "0");
        my $swissprotDBNameUnc = extract_name($configHash->{'swissprotDBLink'}, "gz");
         dl_and_extract($configHash->{'swissprotDBLink'},$swissprotDBNameCom,$swissprotDBNameUnc, $dbDataFolder);
      }
      &run_blast_indexing($configHash->{'swissProtDB'},"blastx",$configHash->{'makeblastdbPath'});
      #Check if the indexing has been good
      if (check_makeblastdb_out($configHash->{'swissProtDB'},$configHash->{'blastdbcmdPath'}) == 1){ 
        #table_in_db_ok(extract_name($configHash->{'swissProtDB'},0));
        print "Swiss-Prot db indexed successfully!\n";
      }else{
        print "There had been problems indexing ".$configHash->{'swissProtDB'}."\n";
      }
    }else { print extract_name($configHash->{'swissProtDB'},0)." already indexed. Going to the next..\n";}
    
    $partDuration = time-$partTime;
    print "Finished the indexing in ".$partDuration." seconds\n\n";
    $partTime = time;
    
    print "Trembl indexed database:\n";
    if ( indexed_db_present(extract_name($configHash->{'tremblDB'},0),\@blastIndexedExts) == 0 ){
      print "Indexing ".$configHash->{'tremblDB'}."\n";
      #If it is defined the link to trembl then download and extract it
      if (defined $configHash->{'tremblDBLink'}){
        my $tremblDBNameCom = extract_name($configHash->{'tremblDBLink'}, "0");
        my $tremblDBNameUnc = extract_name($configHash->{'tremblDBLink'}, "gz");
        dl_and_extract($configHash->{'tremblDBLink'},$tremblDBNameCom,$tremblDBNameUnc, $dbDataFolder); 
      }
      &run_blast_indexing($configHash->{'tremblDB'},"blastx",$configHash->{'makeblastdbPath'});
      #Check if the indexing has been good
      if (check_makeblastdb_out($configHash->{'tremblDB'},$configHash->{'blastdbcmdPath'}) == 1){ 
        #table_in_db_ok(extract_name($configHash->{'tremblDB'},0));
        print "Trembl db indexed successfully!\n";
      }else{
        print "There had been problems indexing ".$configHash->{'tremblDB'}."\n";
      }
    }else { print extract_name($configHash->{'tremblDB'},0)." already indexed. Going to the next..\n";}
    
    $partDuration = time-$partTime;
    print "Finished the indexing in ".$partDuration." seconds\n\n";
    $partTime = time;
  }
  if ( $configHash->{'dbInUse'} eq 'uniref' ){
    print "Swiss-Prot indexed database:\n";
    if ( indexed_db_present(extract_name($configHash->{'swissProtDB'},0), \@blastIndexedExts) == 0 ){
      print "Indexing ".$configHash->{'swissProtDB'}."\n";
      #If it is defined the link to swiss-prot then download and extract it
      if (defined $configHash->{'swissprotDBLink'} ){
        my $swissprotDBNameCom = extract_name($configHash->{'swissprotDBLink'}, "0");
        my $swissprotDBNameUnc = extract_name($configHash->{'swissprotDBLink'}, "gz");
         dl_and_extract($configHash->{'swissprotDBLink'},$swissprotDBNameCom,$swissprotDBNameUnc, $dbDataFolder);
      }
      
      &run_blast_indexing($configHash->{'swissProtDB'},"blastx",$configHash->{'makeblastdbPath'});
      #Check if the indexing has been good
      if (check_makeblastdb_out($configHash->{'swissProtDB'},$configHash->{'blastdbcmdPath'}) == 1){ 
        #table_in_db_ok(extract_name($configHash->{'swissProtDB'},0));
        print "Swiss-Prot db indexed successfully!\n";
      }else{
        print "There had been problems indexing ".$configHash->{'swissProtDB'}."\n";
      }   
    }else { print extract_name($configHash->{'swissProtDB'},0)." already indexed. Going to the next..\n";}
    
    $partDuration = time-$partTime;
    print "Finished the indexing in ".$partDuration." seconds\n\n";
    $partTime = time;
    
    print "Uniref indexed database:\n";
    #The indexed db should be present with all of its files
    if (  indexed_db_present(extract_name($configHash->{'unirefDB'},0), \@blastIndexedExts) == 0 ){
      print "Indexing ".$configHash->{'unirefDB'}."\n";
        #If it is defined the link to trembl then download and extract it
      if (defined $configHash->{'unirefDBLink'}){
        my $unirefDBNameCom = extract_name($configHash->{'unirefDBLink'}, "0");
        my $unirefDBNameUnc = extract_name($configHash->{'unirefDBLink'}, "gz");
        dl_and_extract($configHash->{'unirefDBLink'},$unirefDBNameCom,$unirefDBNameUnc, $dbDataFolder); 
      }
       
      &run_blast_indexing($configHash->{'unirefDB'},"blastx",$configHash->{'makeblastdbPath'});
      #Check if the indexing has been good
      if (check_makeblastdb_out($configHash->{'unirefDB'},$configHash->{'blastdbcmdPath'}) == 1){ 
        #table_in_db_ok(extract_name($configHash->{'unirefDB'},0));
        print "Uniref db indexed successfully!\n";
      }else{
        print "There had been problems indexing ".$configHash->{'unirefDB'}."\n";
      }     
    }else { print extract_name($configHash->{'unirefDB'},0)." already indexed. Going to the next..\n";}
    
    $partDuration = time-$partTime;
    print "Finished the indexing in ".$partDuration." seconds\n\n";
    $partTime = time;
  }
  
 
  print "ncRNA indexed database:\n";
  if ( indexed_db_present(extract_name($configHash->{'rfamDB'},0),\@ncRNAExts) == 0 ){
    print "Indexing ".$configHash->{'rfamDB'}."\n";
    &run_blast_indexing($configHash->{'rfamDB'},"blastn",$configHash->{'makeblastdbPath'});
    #Check if the indexing has been good
    if (check_makeblastdb_out($configHash->{'rfamDB'},$configHash->{'blastdbcmdPath'}) == 1){ 
        #table_in_db_ok(extract_name($configHash->{'rfamDB'},0));
        print "ncRNA db indexed successfully!\n";
    }else{
        print "There had been problems indexing ".$configHash->{'rfamDB'}."\n";
    }     
  }else { print extract_name($configHash->{'rfamDB'},0)." already indexed. Going to the next..\n";}  
  
  $partDuration = time-$partTime;
  print "Finished the indexing in ".$partDuration." seconds\n\n";
  $partTime = time;
  
  if ($configHash->{'keepTempFiles'} ne "YES"){
    delete_directory($configHash->{'tempTableFolder'});
  }
  
	#Calculate the time needed for all the processes
	my $globalDuration = time - $globalStart;
  print "\n(".scalar(localtime).") Annocript has finished the protein database creation in $globalDuration\n";#LOG
} 

=head2 indexed_db_present

 Title   : indexed_db_present
 Usage   : indexed_db_present( -dbToSearch => name of the db that is searching,
                                 -extensions  => an array with the extensions to search,
			       );

 Function: this function checks if a database for the BLAST is consistent

 Returns : nothing

=cut
 sub indexed_db_present {
   my $dbToSearch = shift;
   my $extensions = shift;
   
   #Needed variables
   my $out = '.out';
   my $descrWord = 'Volumes:';
   my $countFiles = 0;  
   my $shouldBe = 0; 
   my $retVal = 0;
   
   #my $dbName = extract_name($dbToSearch,0);
   my $dbName = $dbDataFolder."/".$dbToSearch;
   #go in the dbData folder 
  
   #for the corresponding database we are searching for, open the file with .out extension
   if ( -e $dbName.$out ){
     open(VOL_FILE,"<".$dbName.$out) or die "Cannot open file ".$dbName.$out.". Check you permissions.\n"; 
     my $volReached = 0;
     #print "Reading the file $dbName$out\n";
     #there you should search for the line where is present "Volumes:"
     while (my $line = <VOL_FILE>){
      if ( $line =~ /$descrWord/ ){
        while (my $line1 = <VOL_FILE>){
           #then go to the next lines through the end of the file and read all the path it created
          #for each of the lines (and so the paths) you should check that every extension is there
          #to do it we can call the function indexing_files_present        
          chop($line1);#Remove the \n
          $line1 =~ s/\s//g;
          #print "Checking for: $line1...\n";
          #print "db_creator: The extensions are: ";
          #print_array($extensions);
          #print "\n";
          
          if (indexing_files_present($line1,$extensions) == 1){
             $countFiles +=  scalar(@$extensions);
          }
          $shouldBe += scalar(@$extensions);
        }
      }  
     
     }
     print "Found $countFiles indexing files and they should be: $shouldBe\n";#DEBUGCODE
     close(VOL_FILE);
   #to understand if everything is there the count of files found should be numPaths*numExtensions
   #if it is return true
   if ($countFiles == $shouldBe){
     $retVal = 1; 
     print "\n"
    }else{
      print "There are some problem with $dbName indexing. Annocript will do it again...\n";
    }  
     
   }else{
     print "Cannot open ".$dbName.$out." The indexing of $dbToSearch will start...\n";
   }
  return $retVal;
}

 
=head2 getRfamAndCdd

 Title   : getRfamAndCdd
 Usage   : getRfamAndCdd();

 Function: downloads the CD and RFAM databases, uncompat them and put in the data folder

 Returns : nothing

=cut
 sub getRfamAndCdd {
  my $cdDBLink = $configHash->{'cdDBLink'};
	
  my $cdDBNameCom = extract_name($cdDBLink,"0");
  my $cdDBNameUnc = $configHash->{'cdDB'};
  #Here we obtain the absolute name of the CD database 
  
  my $fileToCheck  = $cdDBNameUnc.".rps";#We check for the .rps file to start
   
  dl_and_extract($cdDBLink,$cdDBNameCom,$fileToCheck, $dbDataFolder); 
  
  #Here we obtain the absolute name of the rFAM database
  my $rfamDBLink = $configHash->{'rfamDBLink'};
  
  my $rfamDBNameCom = extract_name($rfamDBLink,"0");
  my $rfamDBNameUnc = extract_name($rfamDBLink,"gz");
     
  dl_and_extract($rfamDBLink, $rfamDBNameCom, $rfamDBNameUnc, $dbDataFolder);
  
	#Here we merge the ribosomial sequences with these other sequences from a file built by us
	#this file is contained in the Annocript folder
	my $riboDl = $dbDataFolder."/".$rfamDBNameUnc;
	my $ourRiboSeqsNameCom = $configHash->{'usefulFolder'}."/".$configHash->{'ourRiboSeqs'};
	
  #Here we uncompress and get the name of our RNA sequences
  my_extract_file($ourRiboSeqsNameCom,$dbDataFolder);
  my $ourRiboSeqsNameUnc = $dbDataFolder."/".extract_name($ourRiboSeqsNameCom,"gz");
  
  #print "db data folder is: $dbDataFolder and conf->rfam: ".$configHash->{'rfamDB'}."\n";#DEBUGCODE
	join_files_with_cat($riboDl,$ourRiboSeqsNameUnc,$configHash->{'rfamDB'});
	
	#delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$cdDBNameCom);
    delete_file($dbDataFolder."/".$rfamDBNameCom);
    delete_file($riboDl);
    delete_file($ourRiboSeqsNameUnc);#Our uncompressed sequences. The file is in the data folder.
	}
	
}
 

=head2 createGOTable

 Title   : createGOTable
 Usage   : createGOTable( -GOTermsLink => link to a file with GO terms and descriptions,
                          -uniprotWebUser, uniprotWebPass  =>  user-id and password to give to Uniprot to downloadd,
                          -goTableHash => reference to an hash to fill with GO ids from the database
                          - go_create => use this variable to ask to the function if the table should be created or not 
			       );

 Function: creates the GOTable with each of the GO terms, definitions and divisions. It  parses a file with lines
          like this: GO:0000050\tGO:0006594\sGO:0006871\t urea cycle\tP\n


 Returns : nothing

=cut
sub createGOTable{
  #Input Variables
  my $GOTermsLink  = shift;
  my $uniprotWebUser = shift;
  my $uniprotWebPass = shift;
  my ($goTableHash) =shift;
  my $go_create = shift;
  
  
  my $GOtermsFile = extract_name($GOTermsLink,"0");
  my $goTable = 'goTable'; #Name of the goTable
  my $goTablePath =  $configHash->{'tempTableFolder'}."/".$goTable;
    
  #Download file if it doesn't exists
  if (not (-e $dbDataFolder."/".$GOtermsFile) ){
    print "\n Downloading File: ".$GOTermsLink."\n";
    download_file($GOTermsLink,$dbDataFolder,$uniprotWebUser,$uniprotWebPass);
  }else{
    print "File ".$dbDataFolder."/".$GOtermsFile." already downloaded...\n";
  }
    
  open(GOFILE,$dbDataFolder."/".$GOtermsFile);
  if ($go_create == 1){
    open(GOTABLE,">$goTablePath");
  }
     
  #print "...creating a file and an hash with GO information...";
  my $go_id = 1;
  my $goName = "";
  my $definition = "";
  my $division = "";
  my $secPart = "";
  
  #read the GO definitions file
  while(my $row = <GOFILE> ) {
    #The first lines are descriptions and they start with '?'
    if ($row =~ /^!/){
      #DONOTHING
    }
    else{
      
      #$row =~ /^(\S+)\t(.*)?/;
      my @strings = split (/\t/, $row);
      $goName = $strings[0];
      $secPart = $strings[1];
      $definition = $strings[2];
      $division = $strings[3];
      #print "goname: ".$goName; #DEBUGCODE
       if ($secPart){
         #Sometimes in the second part (Separated by tab) other GO terms are there
         #and we should take them all with the same descriptions
         my @keys = split (/\s/, $secPart);
                  
         foreach my $key (@keys){
           #print " - ".$key." - ";#DEBUGCODE
           if ($go_create == 1){
              print GOTABLE join("\t",$go_id, $key, $definition, $division."\n");
           }
           $$goTableHash->{$key} = $go_id;
           #print  join("\t",$go_id, $key, $definition, $division."\n");#DEBUGCODE
           #print "\n";#DEBUGCODE
           $go_id++;
          }
        }
        #One term is always there so... take it!
        if ($go_create == 1){
          print GOTABLE join("\t",$go_id, $goName, $definition, $division."\n");
        }
        $$goTableHash->{$goName} = $go_id;
        #print  join("\t",$go_id, $goName, $definition, $division."\n");#DEBUGCODE
        $go_id++;
        #print "\n";#DEBUGCODE
    }
  }
    
  close(GOFILE);
  if ($go_create == 1){
    close(GOTABLE);
  }
  if ($go_create == 1){
    print " uploading in db...\n";
    table_2_db($goTable);
  }
   
  #delete_file($goTablePath);
 
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$GOtermsFile);
  }
  #rint "....DONE!(GOTABLE Upload completed)\n ";
}


=head2 createPathwayTableAndMapping

 Title   : createPathwayTableAndMapping
 Usage   : createPathwayTableAndMapping( -pathwaysLink => link on thw web to the file with Pathways where to fetch them,
                                 -mappingTableHash  => reference to an hash where the ids for the pathways will be put,
                                 -uniTable => Uniref or Unirefkb, depending on when it is called
                                 -uniprotIdHash => an hash that contains database ids for each Uniprot id
			       );

 Function: downloads the file with the pathways from Uniprot. Then extracts from the lines the pathways and the corresponding
            Uniprot IDs. It creates a table with all the descriptions with an IDs and fills the mappings_hash with the Pathways ids
            corresponding to each Uniprot identifier. It works only for UniPathways


 Returns : nothing

=cut
sub createPathwayTableAndMapping {
  #Input Variables
  my $pathwaysLink  = shift;#Link to the file to download
  my ($mappingTableHash) = shift;#The hash that has to be filled with this new mapping
  my $uniTable = shift;#Uniref or Uniprot_KB database: the Uniprot id will change accordingly
  my $uniprotIdHash = shift;#Hash with uniprot Ids mapped with DB ids
  
  #Extract the name to see if it is already there
  my $fileName = extract_name ($pathwaysLink,"0");
  my $pathwaysTableName = 'pathwaysTable';
  my $pathwaysTablePath =  $configHash->{'tempTableFolder'}."/".$pathwaysTableName;
  
  #Download of the file with Pathways
  if (not (-e $dbDataFolder."/".$fileName) ){
    print "\n Downloading File ".$pathwaysLink."...";
    download_file($pathwaysLink,$dbDataFolder);
    print "...DONE!\n";
  }else{   print "File ".$pathwaysLink." already downloaded...\n"; }
  
  #Here we create a log file where to write informations on 
  #what are the genes that haven't been found in UNIPROT for the GO terms
  print "If Annocript can't find some Uniprot ID will write in ".$dbDataFolder."/".$configHash->{'uniprotGenesNotPresent'}."\n";
  open (UGLOG, ">>".$dbDataFolder."/".$configHash->{'uniprotGenesNotPresent'}) 
  or print "ERROR: Cannot open log file for missing genes in UNIPROT. Maybe you need more permissions.\n";
  
  my $pathwaysFile = $dbDataFolder."/".$fileName;
  
  my $inserted = 0; #inserted elements
  
  
  # Open the pathways file
  open(IN,"<".$pathwaysFile);
  open (PW_TABLE, ">".$pathwaysTablePath);
  
  # Counter for parsed rows
  my $parsedRows = 0;

  # Check the initialization of the href 
  my $ini = 0;

  # Counter for the number of pathways parsed
  my $pathwayId = 1;

  # The row to parse
  my $row;

  # Start the parsing!
  while($row = <IN>) {

    # Row number
    $parsedRows ++;

    # Jump the first 12 lines (TO CHECK EVERY NEW RELEASE!)
    next unless $parsedRows > 12;

    # Jump if no word characters are present in the row
    next unless $row =~ /\w+/;
    chomp($row);

    # This will be the end of the file
    # (TO CHECK EVERY NEW RELEASE!)
    if($row =~ /^Copyrighted/) {
      #print "Found last row:\n$row\n";
      #print "LAST: $id\n";
    }

    # If the row starts with text then we are at the
    # beginning of a new pathway, therefore we will
    # store the collected info for the previous pathway
    # if the initialization has been done
    if($row =~ /^\w+/) {
      if($ini) {
        $pathwayId ++;
        print PW_TABLE "$pathwayId\t$row\n";
      }

      # If the initialization has not been done we are
      # at the beginning therefore this is the first pathway
      else {
        $ini ++;
        print PW_TABLE "$pathwayId\t$row\n";
      }
    }
    
    # If the row does not start with character but with spaces
    # we are in a row containing the ids associated to the
    # current pathway
    else {
      my @sp = split(/\,/,$row);
      foreach my $sp(@sp) {
        if($sp =~ /\((.+)\)/) {
            #Take the uniprot id
            my $unipId = $1;
            
            #If the database is with uniref we must add the 'uniref' word near the id
            if ( $uniTable eq 'uniref' ){
              $unipId = getUnirefType($configHash->{'unirefDBLink'}).'_'.$unipId;
            }
           #Usage of a choice by Annocript to use or not an hash that takes a lot of RAM
           #or to do queries to the database 
            my $unipIndex = '';      
            if ( $configHash->{'useHugeRAM'} eq 'NO'){
              $unipIndex = fetchUniprotId($unipId,$uniTable);  
            }elsif ((defined $uniprotIdHash->{$unipId}) or (exists $uniprotIdHash->{$unipId})){
              $unipIndex =  $uniprotIdHash->{$unipId};
            } 
  
            if ( $unipIndex ne ''){   
                $inserted++; 
                if ( (defined $$mappingTableHash->{$unipIndex}) and ( exists $$mappingTableHash->{$unipIndex}->{'p'}) ){ 
                #Extract all the actual pathways ids and control if the current is present
                #This control is pratically useless because we are flowing the file of pathways and they are sorted
                my @terms = split(",",$$mappingTableHash->{$unipIndex}->{'p'});
                pop(@terms);#This is to erase the last element that is simply a space but the split takes it
                if ( scalar(grep { /\b$pathwayId\b/} @terms) == 0){
                  $$mappingTableHash->{$unipIndex}->{'p'}.=$pathwayId.",";
                  #print "Adding $pathwayId..\n";
                }
                    
              }else{
                $$mappingTableHash->{$unipIndex}->{'p'}.=$pathwayId.",";
                #Don't do any CHOMP. I'll do later when I'll write mapping on the tables
                #print "Adding $pathwayId..\n";
              }  

            }else {	print UGLOG "WARNING (in createPathwayTableAndMapping): there is no $unipId sequence in the uniprotId hash\n";}
        }
        #else {
          #print "ERROR: cannot find SP id for element $sp\n";
        #}
      }
    }
  }
  
  #Close files
  close(IN);
  close(PW_TABLE);
  
  table_2_db($pathwaysTableName);
  #delete_file($pathwaysTablePath);#Tables will be deleted all together
  
  #delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$fileName);
  }
}

=head2 createUniprotKBDB

 Title   : createUniprotKBDB
 Usage   : createUniprotKBDB( -swissprotDBLink => link to the swiss-prot database,
                              -tremblDBLink => 
			       );

 Function: populates a table named UNIPROT with this fields: ( id INTEGER, seqId varchar, dbId INTEGER, descId varchar,
                                        OSId INTEGER, GNId INTEGER, PEId INTEGER, SVId INTEGER)" );
            It downloads both the swiss prot and trembl databases (or exclusively, depending from the link given in input.
            Then uncompress and merge them in a unique file. From this file the headers are extracted and parse with
            some regular expressions.
            Some tables are filled and finally uploaded in the database.
            
            All the file used will be erased. A log file is created with things not found in the fasta database file.

 Returns : nothing

=cut
sub createUniprotKBDB{
  
	#Input variables
  my $swissprotDBLink = shift;
  my $tremblDBLink = shift;
  my $versionDBLink = shift;
  my $uniprotWebUser = shift;
  my $uniprotWebPass = shift;
  my $parsingFilePath = shift;
  
  my $uniprotKBDB = $configHash->{'uniprotKBDB'};  
  $parsingFilePath = $dbDataFolder."/".$parsingFilePath;#Put it in the data folder
  
  #Indexes
  my $idid = 1;
  my $descid = 1;
  my $orgid = 1;
  my $gnid = 1;
  my $peid = 1;
  my $svid = 1;
  my $dbid = 1;
  
  #Temporary variables for fill up the tables
  my $organism = undef;
  my $gn = undef;
  my $pe = undef;
  my $sv = undef;
  my $preid = undef; 
  my $description = undef;

	#File used as table for database (MAYBE IS BETTER IF THEY ARE GLOBAL IN THE TOP OF THE PAGE)
	my $descType = 'uniprotkbDescType';
	my $dbType = 'dbType';
	my $OSType = 'OSType';
	my $GNType = 'GNType';
	my $SVType = 'SVType';
	my $PEType = 'PEType';
	my $uniprotkb = 'uniprotkb';
	
  #Here we create a log file where to write informations on 
  #what are the missing information in the headers of TREMBL and SWISSPROT
  my $logFile = $dbDataFolder."/".$uniprotkb."_".$configHash->{'headParsNAValues'};
  print "A log of the operation will be written in $logFile\n";
  open (HPLOG, ">$logFile") 
  or print "ERROR: Cannot create log file for missing genes in UNIPROT. Maybe you need more permissions.\n";
  
  my $uniprotkbTablePath =  $configHash->{'tempTableFolder'}."/".$uniprotkb;
  #Open the file to write the uniprot table
	open(UNI,">$uniprotkbTablePath");
	
  my $swissprotDBNameCom = undef;
  my $swissprotDBNameUnc = undef;
  
  #If it is defined the link to swiss-prot then download and extract it
  if (defined $swissprotDBLink ){
    $swissprotDBNameCom = extract_name($swissprotDBLink, "0");
    $swissprotDBNameUnc = extract_name($swissprotDBLink, "gz");
     dl_and_extract($swissprotDBLink,$swissprotDBNameCom,$swissprotDBNameUnc, $dbDataFolder);
  }
  
  my $tremblDBNameCom = undef;
  my $tremblDBNameUnc = undef;
 
  #If it is defined the link to trembl then download and extract it
  if (defined $tremblDBLink){
    $tremblDBNameCom = extract_name($tremblDBLink, "0");
    $tremblDBNameUnc = extract_name($tremblDBLink, "gz");
    dl_and_extract($tremblDBLink,$tremblDBNameCom,$tremblDBNameUnc, $dbDataFolder); 
  }

  #If they both are present, we concatenate the two file in one only to obtain an unique database
 if ( (defined $tremblDBLink)  and  (defined $swissprotDBLink) ) {
      print "Creating a file with all the headers in $parsingFilePath (it makes the db construction faster but it will need a while)...";
       my $file =  $dbDataFolder.'/'.$swissprotDBNameUnc;
      try_exec_command( qq(fgrep '>' $file > $parsingFilePath) ) == 1
       or die "ERROR [$?]: an error occurred while extracting headers with fgrep from $file to $parsingFilePath: ?"; 
       $file =   $dbDataFolder.'/'.$tremblDBNameUnc;
       try_exec_command( qq(fgrep '>' $file  >> $parsingFilePath) ) == 1
       or die "ERROR [$?]: an error occurred while appending headers with fgrep from $file  to $parsingFilePath: ?"; 
  }elsif ( !(defined $tremblDBLink)  and  (defined $swissprotDBLink) ){
     $uniprotKBDB = $dbDataFolder."/".$swissprotDBNameUnc;
      #Extracting the headers of the complete file so that it can be parsed more rapidly  
      print "Creating a file with all the headers from $uniprotKBDB (it makes the db construction faster but it will need a while)...";
      try_exec_command( qq(fgrep '>' $uniprotKBDB > $parsingFilePath) ) == 1
       or die "ERROR [$?]: an error occurred while extracting headers with fgrep from $uniprotKBDB to $parsingFilePath: ?";
    }elsif ( (defined $tremblDBLink)  and  !(defined $swissprotDBLink) ){
     $uniprotKBDB = $dbDataFolder."/".$tremblDBNameUnc;
      #Extracting the headers of the complete file so that it can be parsed more rapidly  
      print "Creating a file with all the headers from $uniprotKBDB (it makes the db construction faster but it will need a while)...";
      try_exec_command(  qq(fgrep '>' $uniprotKBDB > $parsingFilePath) ) == 1
       or die "ERROR [$?]: an error occurred while extracting headers with fgrep from $uniprotKBDB to $parsingFilePath: ?";
    }else{
        die "ERROR [?]: both swiss prot and trembl are not defined: ?\n";
      } 
  print "...DONE!\n"; 
  
  print ("Creating UNIPROT table...");
  open(IN,$parsingFilePath);#OPEN PARSE FILE
  
  #Parsing of the file  
  while(my $row = <IN>) {
    if ($row =~ /^\>(\S+)\s(.*)\sOS=/){
      $preid = $1;
      $description = $2;
    }else{ print HPLOG "ERROR in the format of file  while reading $parsingFilePath at line $idid there's no organism!\n";}
      
    if ($row  =~  /\sOS\=(.+?) [A-Z]{2}=(\S+) /){
      $organism=$1;
      #print "OS: ".$organism."\n";
    }  
    elsif($row =~ /\sOS\=(.+?)/){
      $organism = $1;
    }
    else {
      $organism = 'NA';
      print HPLOG "$idid\tOS\n";
      print STDERR "\nERROR: Cannot find organism for $preid\n";
    }
      
      
    if ($row =~ / GN=(\S+)/){ $gn = $1; } else {$gn = "NA";print HPLOG "$idid\tGN\n";}#print "Missing $idid\tGN\n";}#
    #print "GN ".$gn."  \n";
    if ($row =~ / PE=(\S+)/){ $pe = $1; } else {$pe = "NA";print HPLOG "$idid\tPE\n";}#print "Missing $idid\tPE\n";}#DEBUGCODE
    #print "PE ".$pe."  \n";
    if ($row =~ / SV=(\S+)/){ $sv = $1; } else {$sv = "NA";print HPLOG "$idid\tSV\n";}#print "Missing $idid\tSV\n";}#DEBUGCODE
    #print "SV ".$sv."  \n";
       
    #Insert the element just extracted to the hash for each type
    addToHash($descType,$description,\$descid);
    addToHash($OSType,$organism,\$orgid);
    addToHash($GNType,$gn,\$gnid);
    addToHash($PEType,$pe,\$peid); 
    addToHash($SVType,$sv,\$svid) ; 
     
    my @id = split(/\|/,$preid);
    my $seqId = $id[1];
    my $db = $id[0];
    
    addToHash($dbType,$db,\$dbid);
      
    print UNI join("\t",$idid,$seqId,$info->{$dbType}->{$db},$info->{$descType}->{$description},
                     $info->{$OSType}->{$organism},$info->{$GNType}->{$gn},$info->{$PEType}->{$pe},$info->{$SVType}->{$sv})."\n";
    #print  join("\t",$idid,$seqId,$info->{$dbType}->{$db},$info->{$descType}->{$description},
                    #$info->{$OSType}->{$organism},$info->{$GNType}->{$gn},$info->{$PEType}->{$pe},$info->{$SVType}->{$sv})."\n";#DEBUGCODE
     
    $idid++;
    #system "clear";#DEBUGCODE
    #print "\n".scalar(localtime)." Sequence: ".$idid;#DEBUGCODE
  }
  close(UNI);
  close(HPLOG);
    
  #Now write the hashes on file-tables
  print "...DONE! \n Writing the hashes on file-tables....";
  #print Dumper $info;
  
  sortAndWrite($descType);#,$descType);
  print " \n..descriptions table: OK!";
  sortAndWrite($OSType);#, $OSType);
  print " \n....organisms table: OK!";
  sortAndWrite($GNType);#,$GNType);
  print " \n......GN table: OK!";
  sortAndWrite($PEType);#, $PEType);
  print "\n........PE table: OK!";
  sortAndWrite($SVType);#,$SVType);
  print "\n..........SV table OK!";
  sortAndWrite($dbType);#, $dbType);
  print "\n............db types table:OK!\n";
    
  #The table just created can be uploaded in our database using a simple mysql command 
  print "\n Now uploading tables to database: ".$database."...\n";
  
  #Upload tables to db
  table_2_db($OSType);
  table_2_db($descType);
  table_2_db($dbType);
  table_2_db($PEType);
  table_2_db($GNType);
  table_2_db($SVType);
  table_2_db($uniprotkb);
  table_in_db_ok($uniprotkb);
       
  #And finally delete all the tables used that are no more useful
  #delete_file($uniprotkbTablePath);
  #delete_file($configHash->{'tempTableFolder'}."/".$OSType);
  #delete_file($configHash->{'tempTableFolder'}."/".$dbType);
  #delete_file($configHash->{'tempTableFolder'}."/".$PEType);
  #delete_file($configHash->{'tempTableFolder'}."/".$GNType);
  #delete_file($configHash->{'tempTableFolder'}."/".$SVType);
  #delete_file($configHash->{'tempTableFolder'}."/".$descType);
   
  #Delete the log file if it is empty
  if ( -s ($logFile) == 0){
    print "No errors in the parsing. $logFile will be erased.\n";
    delete_file($logFile);
  }
  
  #delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    if (defined ($swissprotDBNameCom) ){delete_file($dbDataFolder."/".$swissprotDBNameCom)};
    if (defined ($tremblDBNameCom) ){delete_file($dbDataFolder."/".$tremblDBNameCom) };
    if (defined ($parsingFilePath) ){delete_file($parsingFilePath)};
  }

}

=head2 createUnirefDB

 Title   : createUnirefDB
 Usage   : createUnirefDB( -db => 'databaseName.db',
			       );

 Function: populates a table named UNIREF with this fields: ( id INTEGER, seqId varchar, n INTEGER, descId varchar,
                                        taxId INTEGER, RepId varchar)" );
            It downloads the uniref database From this file the headers are extracted and parsed with
            some regular expressions.
            Some tables are filled and finally uploaded in the database.
            
            All the files used will be erased. A log file is created with things not found in the fasta database file.
 Returns : nothing

=cut
sub createUnirefDB{
  
	#Input variables
  my $unirefDBLink = shift;
  my $versionDBLink = shift;
  my $uniprotWebUser = shift;
  my $uniprotWebPass = shift;
  my $parsingFilePath = shift;
  
  $parsingFilePath = $dbDataFolder."/".$parsingFilePath;#Put it in the data folder
  
  #Indexes
  my $idid = 1;
  my $descid = 1;
  my $taxid = 1;
  
  #Temporary variables for fill up the tables
  my $taxonomy = undef;
  my $seqId = undef;
  my $description = undef;
  my $n = undef;
  my $repId = undef;

	#File used as table for database
	my $descType = 'unirefDescType';
	my $taxType = 'taxType';
	my $uniref = 'uniref';
	
  my $unirefTablePath =  $configHash->{'tempTableFolder'}."/".$uniref;
  #Open the file to write the uniprot table
	open(UNI,">$unirefTablePath");
	
  #print "Creating UNIPROT db...\n";
   
  #Here we create a log file where to write informations on 
  #what are the missing information in the headers of UNIREF
  my $logFile = $dbDataFolder."/".$uniref."_".$configHash->{'headParsNAValues'};
  print "A log of the operation will be written in $logFile\n";
  open (HPLOG, ">$logFile") 
  or print "ERROR: Cannot create log file for missing genes in UNIREF. Maybe you need more permissions.\n";
  
  my $unirefDBNameCom;
  my $unirefDBNameUnc;
    #If it is defined the link to trembl then download and extract it
  if (defined $unirefDBLink){
    $unirefDBNameCom = extract_name($unirefDBLink, "0");
    $unirefDBNameUnc = extract_name($unirefDBLink, "gz");
    dl_and_extract($unirefDBLink,$unirefDBNameCom,$unirefDBNameUnc, $dbDataFolder); 
  }
  
  
  print "Creating a file with all the information about the sequences. (This file is needed to have a faster process)....PLEASE WAIT...\n";
  
  $unirefDBNameUnc = $dbDataFolder."/".$unirefDBNameUnc;
  
  try_exec_command(qq(fgrep '>' $unirefDBNameUnc > $parsingFilePath) ) == 1
   or die "ERROR [$?]: an error occurred while extracting headers with fgrep from $unirefDBNameUnc to $parsingFilePath: ?\n"; 
  
  print "...DONE!\n";
   
  print ("Creating UNIREF table.....please wait\n");
  open(IN,$parsingFilePath);#OPEN PARSE FILE
  
  #Parsing of the header  
  while(my $row = <IN>) {
    if ($row =~ /^\>(\S+)\s(.*)\sn=(\d+)/){
      $seqId = $1;
      $description = $2;
      $n = $3;
    }else{print HPLOG "There is an error on line $idid. Impossible to read safely the description and ID.\n";}
    
    if ($row  =~  /\s[Tt]ax\=(.+?) RepID=(\S+)/){
      $taxonomy = $1;
    } 
    elsif($row =~ /\s[Tt]ax\=(.*)/){
      $taxonomy = $1;
    }
    else {
      $taxonomy = "NA";
      print STDERR "\nERROR: Cannot find tax\n";
      print HPLOG "$idid\ttax\n";
    }
    
    if ($row =~ / RepID=(\S+)/){ $repId = $1; } else {$repId = "NA";print HPLOG "$idid\t repID\n";}
       
    #Insert the element just extracted to the hash for each type
    addToHash($descType,$description,\$descid);
    addToHash($taxType,$taxonomy,\$taxid); 
     
    print UNI join("\t",$idid,$seqId,$info->{$descType}->{$description}, $n, $info->{$taxType}->{$taxonomy},$repId)."\n";
    #print  join("\t",$idid,$seqId,$info->{$descType}->{$description}, $n, $info->{$taxType}->{$taxonomy},$repId)."\n";#DEBUGCODE
     
    $idid++;
    #system "clear";#DEBUGCODE
    #print "\n".scalar(localtime)." Sequence: ".$idid;#DEBUGCODE
  }
  close(UNI);
  close(HPLOG);
    
  #Now write the hashes on file-tables
  print "\n $dbDownloaded table completed! \n Writing the hashes on file-tables....";
  #print Dumper $info;
  
  sortAndWrite($descType);#,$descType);
  print " \n....descriptions table..ok...";
  sortAndWrite($taxType);#, $taxType);
  print " ....taxType table..ok...";
   
  #The table just created can be uploaded in our database using a simple mysql command 
  print "...hash completed!  \n Uploading tables to database: ".$database."...\n";
  
  #Upload tables to db
  table_2_db($taxType);
  table_2_db($descType);
  table_2_db($uniref);
  table_in_db_ok($uniref);
       
  #And finally delete all the tables used that are no more useful
  #delete_file($unirefTablePath);
  #delete_file($configHash->{'tempTableFolder'}."/".$taxType);
  #delete_file($configHash->{'tempTableFolder'}."/".$descType);
   
  #Delete the log file if it is empty
  if ( -s ($logFile) == 0){
    print "No errors in the parsing. $logFile will be erased.\n";
    delete_file($logFile);
  }
  
  #delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$unirefDBNameCom);
    delete_file($parsingFilePath);
  }
  print "..completed!\n";
}


=head2 table_2_db

 Title   : table_2_db
 Usage   : table_2_db( -table => 'table to import'
			       );
             

 Function: 	Uses DBI to import a CSV file in a table
 Returns : nothing

=cut
sub table_2_db{
  my $table = shift;
  
  my $path_to_table = $tempTableFolder."/".$table;#Creating the path to the table
  my $query = "LOAD DATA LOCAL INFILE '".$path_to_table."' INTO TABLE ".$table;
  
  #query_to_db("LOAD DATA INFILE '$table' INTO TABLE $table FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n'");
  unless (query_to_db($query) == 1){
    print "Cannot execute query: $query\n";
    print "Annocript will try to execute with a system call...\n";
    table_2_db_sys($table);
  }
}

      
=head2 table_2_db_sys

 Title   : table_2_db_sys
 Usage   : table_2_db_sys( -table => 'table to import'
			       );
             

 Function: 		Imports a CSV file in a specified table of the database with a system call. It calls the method try_exec
 that will try a set of times.
 
 Returns : nothing

=cut
sub table_2_db_sys{
  my $table = shift;
  
  my $maxTent = 5;
  my $tentatives = 0;
  my $path_to_table = $tempTableFolder."/".$table;#Creating the path to the table
  my $threads = $configHash->{'threads4Parallel'};
  
  
  if ( try_exec_command ("mysqlimport -u$mySqlUser -p$mySqlPass $database $path_to_table --use-threads=$threads --local --fields-terminated-by='\\t' --lines-terminated-by='\\n'" ) < 1){
    #if ( try_exec_command ("mysqlimport -u$mySqlUser -p$mySqlPass $database $table --local " ) < 1){
     #destroy_db($database);
    die "ERROR: an error occurred while importing $table in $database. $? \n";
  }
}

=head2 is_table_in_db

 Title   : is_table_in_db
 Usage   : is_table_in_db( -table => 'table to check'
			       );
             

 Function: Checks if the table "table" is in the database
 
 Returns : 1 if it is there, -1 otherwise

=cut
sub is_table_in_db{
  my $table = shift;
  
  my $retVal = -1;
  
  #Opening database if needed 
	my $dbh;
	if (defined $database){
		# PERL DBI CONNECT AND CREATE TABLE
		$dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
	}
  #Here we open an access to the database to collect informations belonging to the HSPs
  #print "seqID: $seqId\n";
  my @res = ();
  my $select_handle;
  
  #Extracts information from UNIPROT only if $database is defined. (i.e. blastx)	
  if (defined $database){
    #fetch informations from the database
    $select_handle = $dbh->prepare("SELECT * FROM goodTables WHERE tableName='".$table."'");
    die "ERROR [$?]: can't prepare queries; aborting: ?\n" unless defined $select_handle;
    $select_handle->execute() or die $select_handle->errstr;
    @res = $select_handle->fetchrow_array;
    if (scalar(@res) == 0){
      print "There is no $table in $database...creating it...\n";
    }else{
        $retVal = 1;
      }
	 }
   $select_handle->finish;
   	$dbh->disconnect(); 
   return $retVal;
}

=head2 table_in_db_ok

 Title   : table_in_db_ok
 Usage   : table_in_db_ok( -table => 'table to import'
			       );
             

 Function: The table "table" is in the database. We put its name inside 'goodTables' table
 Returns : nothing

=cut
sub table_in_db_ok{
  my $table = shift;
  
  #Do the insert query
  query_to_db("INSERT INTO goodTables (tableName) VALUES ('".$table."')");
  
}


=head2 query_to_db

 Title   : query_to_db
 Usage   : query_to_db( -query => 'the query for the insert'
			       );
             

 Function: The query command is applied
 Returns : nothing

=cut
sub query_to_db{
  my $query = shift;

  #Opening database if needed 
	my $dbh;
	my $dsnToLoadSafely = "dbi:$platform:$database:$host:$port;mysql_local_infile=1";
  my $retVal = 1;
  # PERL DBI CONNECT AND CREATE TABLE
	$dbh = DBI->connect($dsnToLoadSafely, $mySqlUser, $mySqlPass, { RaiseError => 0, AutoCommit => 0 }) 
    or die $DBI::errstr;

  my $select_handle;
  #print "Executing query: $query\n";
  
  #prepares query to the database
  $select_handle = $dbh->prepare ($query);
  die "ERROR [$?]: can't prepare queries; aborting: ?\n" unless defined $select_handle;
  unless (defined $select_handle->execute()){
    $retVal = 0;
  }
   
  $dbh->commit();
	$dbh->disconnect(); 
  return $retVal;
}


=head2 extractEnzimeIdsFile

 Title   : extractEnzimeIdsFile
 Usage   : extractEnzimeIdsFile( -db => 'databaseName.db',
                                -uniGoPath => '/folder/GODBLink'
			       );

 Function: 		Extract enzyme ids from the file with all the uniprot ids. An example of part of the
            enzyme.dat file is
            ID   1.1.1.1
            DE   Alcohol dehydrogenase.
            AN   Aldehyde reductase.
            CA   (1) An alcohol + NAD(+) = an aldehyde or ketone + NADH.
            CF   Zinc or iron.
            CC   -!- Acts on primary or secondary alcohols or hemi-acetalss with very
            PR   PROSITE; PDOC00060;
            DR   P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU;  Q5RBP7, ADH1A_PONAB;
          
            Where we need exclusively the ID, DE (description), DR (uniprot ids related).
            To extract this information for a table we parsed this file.
            
 Returns : Returns a file with 3 columns: enzymeId, Description, UniprotIds corresponding

=cut     
sub extractEnzimeIdsFile{
  my $file = shift;
  
  my $id = undef;
	my $de = undef;
	my $an = undef;
	
  my @strings = ();
	my $uniprotIds = "";
	   
	my $firstId = 0;
	my $idNum = 0;
	my $descNum = 0;
	my $drNum = 0;
	   
	print "Extracting informations about enzymes...\n";
  
  #Creating a path to a file that will contain the table with uniprot ids and enzymes
  my $uniprotEnzymeFilePath = $configHash->{'tempTableFolder'}."/tempEnzymeFile";
	
  open(FILE,$file);
	open (ENZYMEFILE, ">$uniprotEnzymeFilePath");
	  
	while (my $row = <FILE>) {
    #If you find a row with ID, it is the identifier
  	if ($row =~ /^ID\s+(.*)/){
			if($firstId > 0){
				#print "Previous Values: ID ".$id." and descs: ".$de."\n";
				#print "IDS: ".$uniprotIds;
				#print "\n";
				print ENZYMEFILE join("\t",$id,$de,$uniprotIds."\n");
			}
		  $de = "";
			$id = $1;
			$uniprotIds = "";
			$firstId=1;
    }
    
    #The row with the description is with DE
		elsif ($row =~ /^DE\s+(.*)/){
      if ($de eq ''){
        $de = $1;
      }else{
        $de = $de."  ".$1;
      }
					
		$de =~ s/[^A-Za-z0-9\-\(\)\s\_\/\.\+]/ /g;
    
    #The line with DR is the line with the uniprot ids but they must be separated
    }elsif ($row =~ /^DR\s+(.*)/){
			my $line = $1;
			$line =~ s/\s+//;
      #Separation is obtained by taking elements separated by ';'
			@strings  = split (/;/, $line);
				foreach my $string (@strings){
					my @ids = split (/,/ , $string);
					$uniprotIds= $uniprotIds.$ids[0].",";
				}
			}
												     
  }
  
	#print "Previous Values: ID ".$id." and descs: ".$de."\n";
	#print "IDS: ".$uniprotIds;
	#print "\n";
	print ENZYMEFILE join("\t",$id,$de,$uniprotIds);
	#print "..DONE! (A file called $uniprotEnzymeFilePath has been created)\n";
	close(ENZYMEFILE);
	close(FILE);
	return $uniprotEnzymeFilePath;		
}   
 
	

=head2 createCDTable

 Title   : createCDTable
 Usage   : createCDTable( -db => 'databaseName.db',
                           -GOTermsPath  => '/folder/GOTermsFile',
			       );

 Function: creates the CD (conserved domains) Table with each of the domains descriptions. Descriptions
          are taken from a downloaded file. They are 5 columns separated by TAB.
          An hash given in input (cdTableHash) is built using the informations about the domains and
          it will be used to add the informations to the final tables.
 Returns : nothing

=cut
sub createCDTable{
  #Input Variables
  my $tableFile = shift;  
  my $cdTableLink  = shift;
  my $uniprotWebUser = shift;
  my $uniprotWebPass = shift;

  my $compressed = extract_name($cdTableLink,'0');
  my $uncompressed= extract_name($cdTableLink,'gz');#Name for uncompressed file
  
  dl_and_extract($cdTableLink,$compressed,$uncompressed,$dbDataFolder,$uniprotWebUser,$uniprotWebPass);
       
   my $cdTablePath =  $configHash->{'tempTableFolder'}."/".$tableFile;
  open(CDFILE,$dbDataFolder."/".$uncompressed);
  open(CDTABLE,">$cdTablePath");
     
  print "...creating a file with table...";
  my $cd_id = 1;
  my $genericId = "";
  my $cdAccession = "";
  my $cdDescription = "";
  my $pssmLength = "";

  while ( my $row = <CDFILE> ) {
    my @strings = split (/\t/, $row);
    my $genericId = $strings[0];#$1;
    my $cdName = $strings[1];#$2;
    my $cdShortName = $strings[2];
    my $cdDescription = $strings[3];
    my $pssmLength = $strings[4];
       
    print CDTABLE join("\t",$cd_id, $genericId, $cdName, $cdShortName, $cdDescription, $pssmLength."\n");
    $cd_id++;
  }     

  print "..DONE! \n Now uploading in db...";
  
  #imports table to db
  table_2_db($tableFile);
  table_in_db_ok($tableFile);
  
  #Deletes the file-table
  #delete_file($cdTablePath);
   
  #Deletes the downloaded file
  if ($configHash->{'eraseFiles'} eq "YES"){
    if (unlink($dbDataFolder."/".$uncompressed) == 0) {
      print "File deleted successfully.";
    } else {print "File $uncompressed was not deleted.";}
  }
  print "....DONE!(CDTABLE Upload completed)\n ";
  
  close(CDFILE);
  close(CDTABLE);
}


=head2 extractUniprotIdHash

 Title   : extractUniprotIdHash
 Usage   : extractUniprotIdHash( -uniprotIdHash => reference to an hash to fill,
                           -uniTable  => the Uniprot database currently used: uniprotkb or uniref,
			       );

 Function: fills an hash with UniprotId -> DBId by using a set of queries
 Returns : nothing

=cut
sub extractUniprotIdHash{
  my ($uniprotIdHash) = shift;
  my $uniTable = shift;
  
  # CREATING UNIPROT IDs HASH
  # PERL DBI CONNECT AND CREATE TABLE
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
	if ($dbh->err()) { die "$DBI::errstr\n"; }
  
	#my $uniprotIdHash;
	#print "Creating an hash with $uniTable ids...";
  
  #print 'SELECT id, seqId FROM '.$uniTable.';';
  my $select_handle = $dbh->prepare('SELECT id, seqId FROM '.$uniTable.';');
  die "ERROR [$?]: can't prepare queries; aborting: ?" unless defined $select_handle;
  $select_handle->execute();
  while (my $res = $select_handle->fetchrow_hashref){
    $$uniprotIdHash->{$res->{seqId}} = $res->{id};
   # print "res: ". $res->{seqId}."  = ".$res->{id}."\n";
  }
  
  #print Dumper $uniprotIdHash;#DEBUGCODE
  
  $dbh->commit();
  $dbh->disconnect();
}

=head2 checkIdMappingFile

 Title   : checkIdMappingFile
 Usage   : checkIdMappingFile(
			       );

 Function: checks if the column registered in Annocript is equal to the one that Uniprot gave in the mapping file. 

 Returns : 1 if it is ok, 0 otherwise

=cut
sub checkIdMappingFile{
  my $unirefIdPos4Cut = shift;
  my $unirefType = shift;
  my $GODBPath = shift;
  
  my $retVal = 1;
  
  open(GODBFILE,"<$GODBPath");
  my $firstLine = <GODBFILE>;
  if ( (split(/\t/,$firstLine))[$unirefIdPos4Cut] !~ /$unirefType/){
      $retVal = 0;
  }else{
      $retVal = 1;
    }
  close(GODBFILE);
  
  return $retVal;
   
}

=head2 extractSelectedUnirefColumnFromMapping

 Title   : extractSelectedUnirefColumnFromMapping
 Usage   : extractSelectedUnirefColumnFromMapping(
			       );

 Function: extract the column corresponding to the Uniref database in use from the idmapping file 

 Returns : the column searched

=cut
sub extractSelectedUnirefColumnFromMapping{
  my $unirefType = shift;
  my $GODBPath = shift;
  
  my $unirefIdPos4Cut = -1;
  
  open(GODBFILE,"<$GODBPath");
  
  #This loop continues to read lines until the Uniref specific id column is found
  while ( $unirefIdPos4Cut < 0) {
    my $firstLine = <GODBFILE>;
    my @cols = split (/\t/,$firstLine);#Separate per tab
    my $ind= 0;
    for ($ind = 0; $ind < scalar(@cols); $ind++){#Loop across the columns
        if($cols[$ind] =~ /$unirefType/){
            $unirefIdPos4Cut = $ind+1;#Plus 1 because the split starts from 0
        }
    }
  }
    
  close(GODBFILE);
  
  return $unirefIdPos4Cut;
   
}


=head2 extractFilteredMappingFile

 Title   : extractFilteredMappingFile
 Usage   : extractFilteredMappingFile(
			       );

 Function: downloads the file with mapping of Uniprot ids with GO ids and creates a reduced file with only the columns 
            needed. 
            BE CAREFUL!!:
            The columns indexes are hard coded. If Uniprot changes the file, they should be changed accordingly

 Returns : the path to a reduced file with the mapping

=cut
 sub extractFilteredMappingFile{
  #my $GODBLink = shift;
  
  #Positions to use for the CUT program. If the file mapping changes, they should be changed
  my $uniprotIdPos4Cut = 1;
  my $goPos4Cut = 7;
  my $unirefIdPos4Cut;
  
  #Here we get the type of Uniref that we are using to know what column should be cut
  #print "Searching from ".$configHash->{'unirefDBLink'}."\n";#DEBUGCODE
  
  my $unirefType = getUnirefType($configHash->{'unirefDBLink'});
  #print "obtained: $unirefType\n";#DEBUGCODE
  if ( $unirefType eq 'UniRef100'){
      $unirefIdPos4Cut = $configHash->{'UniRef100PosIdMap'};
  }
  if ( $unirefType  eq 'UniRef90'){
      $unirefIdPos4Cut = $configHash->{'UniRef90PosIdMap'};
  }
  if ( $unirefType  eq 'UniRef50'){
      $unirefIdPos4Cut = $configHash->{'UniRef50PosIdMap'};
  }
 
  
  #Extracts the compressed and uncompressed names from the link
  my $GODBLink = $configHash->{'GODBLink'};
  my $GODBNameCom = extract_name($GODBLink,"0");
  my $GODBNameUnc = extract_name($GODBLink,"gz");
   
  
   
  #Extracts the GO table from web
  dl_and_extract($GODBLink,$GODBNameCom,$GODBNameUnc, $dbDataFolder,$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'}); 
  
  #Search the column where it is written the specific Uniref id associated
  $unirefIdPos4Cut = extractSelectedUnirefColumnFromMapping($unirefType,$dbDataFolder."/".$GODBNameUnc);
  
  my $filtMappingFile  = $configHash->{'filtMappingFile'};
  
  
  #if ( checkIdMappingFile($unirefIdPos4Cut,$unirefType,$dbDataFolder."/".$GODBNameUnc) == 1){ 
   if( $unirefIdPos4Cut > 0){
    #With the cut command the columns start from 1
    my $command = "cut -f".$uniprotIdPos4Cut.",".$goPos4Cut.",".$unirefIdPos4Cut." ".$dbDataFolder."/".$GODBNameUnc." > ".$dbDataFolder."/".$filtMappingFile;
    print "Creating a filtered mapping file with the command: $command...\n";
    if ( try_exec_command($command) < 1){
      delete_file($dbDataFolder."/".$filtMappingFile);
      die "Cannot execute $command: ?\n";
    }
  #}else{print "The column $unirefIdPos4Cut in ".$dbDataFolder."/".$GODBNameUnc." does not correspond anymore to $unirefType.".
  #"You can modify the parameters UniRef[VER]PosIdMap in the config_annocript.txt file accordingly or you can contact the creator of Annocript.\n";}
 }else{print "The column $unirefIdPos4Cut in ".$dbDataFolder."/".$GODBNameUnc." does not correspond anymore to $unirefType. Please contact Annocritp creators.\n"} 
 
 
  #delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$GODBNameCom);
    delete_file($dbDataFolder."/".$GODBNameUnc);
  }
  
  return $filtMappingFile;
}

=head2 createMappingTable

 Title   : createMappingTable
 Usage   : createMappingTable( -uniTable => 'name of the uniprot table where we take the ids (uniref, uniprotkb)
                                                his parameter gives the name of the uniprot genes principal table 
                                                but it also can give indication of what we are using',
                                -mappingsTable => 'name of the mappings table'
                                -go_enz_create => with this flag is decided if the GO and ENZYME tables have to be uploaded in db or not
			       );

 Function: until now it creates mapping table for GO terms and Enzyme ones. Functions to call some other mappings
			have to be added here. 
			For GOTerms it is first created the table by downloading the file with mappings. Then parsed and created
			an hash.
			Then an hash of uniprot ids is created i. The presence of this latter can occupy a lot of memory so the function will use queries
      when the memory is less than a certain threshold.
			Then the function to insert GO elements is called with these two hashes.
			Then the function to extract enzyme ids from the file and fill up the mappings table is called.
			These functions can be called separately and exclusively since the inserting system is based 
			on the INSERT-IF-EXISTS type of query.

 Returns : nothing

=cut
sub createMappingTable{
  my $uniTable = shift;#used for both a table and an indication of db used
  my $mappingsTable = shift;
  my $go_create = shift;
  
  #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
	
  #An hash to store all the mapping table
  my $mappingTableHash;
  
  my $mappingsTablePath =  $configHash->{'tempTableFolder'}."/".$mappingsTable;
  print "\n**MAPPING FOR $uniTable**\n";
  
	#CREATING GO TABLE AND ITS HASH 
  my $goTableHash; 
  print "Creating an hash with GO ids from $uniTable...\n";
  createGOTable($configHash->{'GOTermsLink'},$configHash->{'uniprotWebUser'},$configHash->{'uniprotWebPass'},\$goTableHash,$go_create);#Create the Database with all the GO definitions and divisions 
  #print "\n....DONE!\n";
  $partDuration = time-$partTime;
  print "Finished Creating GO table in ".$partDuration." seconds\n\n";
  $partTime = time;
	
	#print Dumper $goTableHash;#DEBUGCODE
  
  #Extracting a reduced file with the GO mapping for the 
  my $filtMappingFile = extractFilteredMappingFile();
  
  my $uniprotIdHash;
  if ( $configHash->{'useHugeRAM'} eq 'YES'){
    extractUniprotIdHash(\$uniprotIdHash,$uniTable);
  }
  
  #print Dumper $uniprotIdHash;#DEBUGCODE

  #GO MAPPING
  print "Reading the GO mapping file and storing terms for $uniTable ids...\n";
 
  if ( $configHash->{'useHugeRAM'} eq 'YES'){
    createUniprotGO($filtMappingFile, $goTableHash,$uniTable,$go_create,\$mappingTableHash, $uniprotIdHash); #Creates an hash table with all uniprot Id and sequence codes
  }else{
    createUniprotGO($filtMappingFile, $goTableHash,$uniTable,$go_create,\$mappingTableHash); #Creates an hash table with all uniprot Id and sequence codes
  } 
  $partDuration = time-$partTime;
  print "...DONE! ..in ".$partDuration." seconds\n";
  $partTime = time;
  
  if ( $uniTable ne 'uniref'){
    print "Reading the ENZYME mapping file and storing terms for $uniTable ids...\n";
    if ( $configHash->{'useHugeRAM'} eq 'YES'){
      createEnzymeTable($configHash->{'enzymeDBLink'},\$mappingTableHash,$uniTable, $uniprotIdHash);
    }else{
      createEnzymeTable($configHash->{'enzymeDBLink'},\$mappingTableHash,$uniTable);
    }
    $partDuration = time-$partTime;
    print "...DONE! ..in ".$partDuration." seconds\n\n";
    $partTime = time;
    
    print "Reading the PATHWAYS mapping file and storing terms for $uniTable ids...\n";
    if ( $configHash->{'useHugeRAM'} eq 'YES'){
      createPathwayTableAndMapping($configHash->{'pathwaysTableLink'},\$mappingTableHash,$uniTable, $uniprotIdHash);
    }else{
      createPathwayTableAndMapping($configHash->{'pathwaysTableLink'},\$mappingTableHash,$uniTable);
    }
    $partDuration = time-$partTime;
    print "...DONE! ..in ".$partDuration." seconds\n\n";
    $partTime = time;
  }
  
  #Here we can start to deallocate these to hashes which are no more useful
  #since everything is contained in the mappingTableHash
  deallocate_hash(\$goTableHash);
  if ( $configHash->{'useHugeRAM'} eq 'YES'){
    deallocate_hash(\$uniprotIdHash);#This is the hughest
  }
  
  ##print Dumper $mappingTableHash;#DEBUGCODE
 
  #Here, after the hash has been filled with mapping informations, we write a table-file to upload in db
  print "Creating a table to be loaded in the db with all information stored...\n";
  
  #For each uniprot ID we will create a line with the corresponding mappings added in the hash
  open(MAPPING_FILE, ">$mappingsTablePath");
  #We will write a line in the table only if something has been added
  my $somethingAdded = 0;
  
  #while (my ($uniprotKey, $value) = each %$mappingTableHash) {  #Medium speed
  #Sorts the hash. This gives speed to the LOAD in db
  #foreach  my $uniprotKey (sort {$mappingTableHash->{$a} <=> $mappingTableHash->{$b}}#Slowest
              #keys %{$mappingTableHash}) {
  foreach my $uniprotKey (keys %$mappingTableHash){#Highest speed
    #The fields are initialized with the SQL null value \N
    my $goTerms = '\N';
    my $enzymeIds = '\N';
    my $pwIds = '\N';
    
    #Initializing an array for the fields to put in the table
    my @fieldsForMapTable = ();
    
    push(@fieldsForMapTable,$uniprotKey);#Adding the key to the hash
    #Add the go terms
    if ( exists $mappingTableHash->{$uniprotKey}->{'g'} ){
      $goTerms = $mappingTableHash->{$uniprotKey}->{'g'};
    	chop($goTerms); 
      $somethingAdded++;
    }
    push(@fieldsForMapTable, $goTerms);

   #Add elements only for swiss-prot
    if ( $uniTable ne 'uniref'){
      #Adding enzyme ids
      if ( exists $mappingTableHash->{$uniprotKey}->{'e'} and $mappingTableHash->{$uniprotKey}->{'e'} ne ' '){
        $enzymeIds = $mappingTableHash->{$uniprotKey}->{'e'};
        chop($enzymeIds);
        $somethingAdded++;
      }
       push(@fieldsForMapTable, $enzymeIds);

      #Adding pathways ids
      if ( exists $mappingTableHash->{$uniprotKey}->{'p'} and $mappingTableHash->{$uniprotKey}->{'p'} ne ' '){
        $pwIds = $mappingTableHash->{$uniprotKey}->{'p'};
        chop($pwIds);
        $somethingAdded++;
      }     
      push(@fieldsForMapTable, $pwIds);
    }
    if ($somethingAdded > 0){
      #push(@fieldsForMapTable,"\n");
      print MAPPING_FILE join("\t",@fieldsForMapTable);
      print MAPPING_FILE "\n";
    }
    $somethingAdded = 0;
  }
  close(MAPPING_FILE);
  
  #deallocate the huge hash
  deallocate_hash(\$mappingTableHash);
  
  print "Loading in DB...\n";
  table_2_db($mappingsTable);#table_2_db takes in input the table name
  #delete_file($mappingsTablePath);#delete_file takes all the path of the file to delete
  
  sleep 5;#Sleeps 5  seconds to be sure the previous process finishes to load the table
  table_in_db_ok($mappingsTable);
    
  #Computing time  
  #print "\n....DONE!\n";
  $partDuration = time - $partTime;
  print "Finished generating and loading $uniTable mapping table in ".$partDuration." seconds\n\n";
  $partTime = time;
}	

=head2 createUniprotGO

 Title   : createUniprotGO
 Usage   : createUniprotGO( -filtMappingFile => 'link to a reduced file with uniprot/uniref ids and corresponding GO terms',
                                -goTableHash => 'hash with the go database ids',
                                -dbInUse => 'uniprot or uniref',
                                -go_create => 'says if the GO table will be created or not',
                                -($mappingTableHash) => 'the hash to fill up for have a mapping',
                                -$uniprotIdHash => 'a reference to an hash with the uniprot/uniref db ids (if used)'
			       );

 Function: fills the mapping table hash with GO terms corresponding to Uniprot/Uniref ids.
					It doesn't insert any NA value and inserts all GOID separated by commas in one line.

 Returns : nothing

=cut
 sub createUniprotGO{
  
  my $filtMappingFile = shift;
  my $goTableHash = shift;
  my $dbInUse = shift;
  my $go_create = shift;
  my ($mappingTableHash) = shift;
  my $uniprotIdHash = shift;
   
    
  my $uni_name;#Index to find the uniprot name if it is Uniref or Uniprot_kb
  my $geneName;
  my $goPosition = 1;#Position in the filtered file for the mapping
  
  #Positions to use for the CUT program. If the file mapping changes, they should be changed
  #my $uniprotIdPos4Cut = 1;
  #my $goPos4Cut = 7;
  #my $unirefIdPos4Cut = 10;
  
  #Assign positions using the filtered file for the mapping
  if($dbInUse eq 'uniprotkb' ){
      $uni_name = 0;
  }elsif($dbInUse eq 'uniref' ){
      $uni_name = 2;
    }else{
      die "ERROR [$?]: the variable $dbInUse has not been defined. Some error is in the code: ?\n";
    }
  
  #Here we create a log file where to write informations on 
  #what are the genes that haven't been found in UNIPROT for the GO terms
  print "A log of the operation will be written in ".$dbDataFolder."/$dbInUse".$configHash->{'uniprotGenesNotPresent'}."\n";
  open (UGLOG, ">".$dbDataFolder."/$dbInUse".$configHash->{'uniprotGenesNotPresent'}) 
  or print "ERROR: Cannot create log file for missing genes in database. Maybe you need more permissions.\n";
  
  my $goId;
  #my $hash_ref;#USELESS
  my $geneNamesNotFound = 0;#Gene names that are not found in the database
 
  #To keep in count the time occurred
  my $partTime = time;
  my $partDuration = undef;
 	
  print "Opening : ".$dbDataFolder."/".$filtMappingFile."\n";
  open(UNIGOFILE,"<".$dbDataFolder."/".$filtMappingFile); 
      
  #print "Updating the mapping hash with GO ids...\n";
  my $rowNum = 0;
  while (my $row = <UNIGOFILE>) {
    #print $rowNum." ";
    $rowNum++;
    #Separate elements by tab
    my @strings = split (/\t/, $row);
    $geneName = $strings[$uni_name]; # The gene name is the first element
    #print "Gene name: $geneName - ";#DEBUGCODE
    
    #If the word contains the \n, remove it. This happens when reading the last field of the file  
    if ($geneName =~ /\n/){
       chop($geneName);
    } 
    #take the string with the goTerms: A POSITION MUST BE chose!
    my $goTermsAss = $strings[$goPosition];
    $goTermsAss =~ s/\s//g;#Erase spaces that are in the string
      
    #Separate GO terms
    my @goTerms = split(";",$goTermsAss);
     
    #If the gene name is found and defined in the hash mapping from the DB         
    if ( ($geneName ne '')){ #and (defined $uniprotIdHash->{$geneName})){
      my $unipIndex = '';
      if ( $configHash->{'useHugeRAM'} eq 'NO'){
        $unipIndex  = fetchUniprotId($geneName,$dbInUse);
      }elsif (defined $uniprotIdHash->{$geneName}){
        $unipIndex = $uniprotIdHash->{$geneName};
      }
      
      if ( $unipIndex ne ''){
        my $goIds = "";
        #For each term put it in the hash on the corresponding uniprot Id
        foreach my $goTerm ( @goTerms){
          if (defined $goTableHash->{$goTerm}){
            #build the go ids sequence
            my $goTableId = $goTableHash->{$goTerm};
              
            #If the corresponding uniprot ID in the hash has not been created yet then we write it 
            if ( (defined $$mappingTableHash->{$unipIndex}) and (exists $$mappingTableHash->{$unipIndex}->{'g'}) ){
                #We build an array with all the terms already inserted
                my @terms = split(",",$$mappingTableHash->{$unipIndex}->{'g'});
                
                #if this term is already there, then do not insert
                if ( scalar(grep {/\b$goTableId\b/} @terms) == 0){
                  $$mappingTableHash->{$unipIndex}->{'g'}.=$goTableId.",";
                }
            }else{
                  $$mappingTableHash->{$unipIndex}->{'g'}.=$goTableId.",";
            }
          }#else {print "ERROR there is no go term $goTerm in the go hash. Check it!";last;}
        }        
      }else{
        $geneNamesNotFound++;
        print UGLOG "WARNING (in createUniprotGO): there is no $geneName in the hash created for the $dbDownloaded ids!!\n";
      }
    }
     
    #$partDuration = time - $partTime;
    #print "For one row of GO terms occurs ".$partDuration." seconds\n";
    #$partTime = time;
  }
  if ($geneNamesNotFound > 0){
      print "Gene names not found: $geneNamesNotFound \n";
  }

   
  #print "..DONE!\n";
  close(UNIGOFILE);
  close(UGLOG);
  
  #Delete the log file if it is empty
  if ( -s ($dbDataFolder."/$dbInUse".$configHash->{'uniprotGenesNotPresent'}) == 0){
    print "No errors in the parsing ".$dbDataFolder."/$dbInUse".$configHash->{'uniprotGenesNotPresent'}." will be erased.\n";
    delete_file($dbDataFolder."/$dbInUse".$configHash->{'uniprotGenesNotPresent'});
  }
}	



=head2 fetchUniprotId
 Title   : fetchUniprotId
 Usage   : fetchUniprotId( -db => 'databaseName.db',
                           -GOTermsPath  => '/folder/GOTermsFile',
			       );

 Function: fetches a db ID given the uniprot ID
 Returns : nothing

=cut
sub fetchUniprotId{
  my $seqId = shift;
  my $uniTable = shift;
  
	# CREATING UNIPROT IDs HASH
  # PERL DBI CONNECT AND CREATE TABLE
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
	if ($dbh->err()) { die "$DBI::errstr\n"; }
  
	my $uniprotId = '';
	#print "Creating an hash with $uniTable ids...";
  
  #print 'SELECT id, seqId FROM '.$uniTable.';';
  my $select_handle = $dbh->prepare("SELECT id FROM $uniTable WHERE seqId='$seqId';");
  die "ERROR [$?]: can't prepare queries; aborting: ?" unless defined $select_handle;
  $select_handle->execute() or die $select_handle->errstr;
  my @res = $select_handle->fetchrow_array;
  if (scalar(@res) > 0){
    $uniprotId = $res[0];
  }#else{ print "WARNING: There is no result for the query (SELECT id FROM $uniTable WHERE seqId='$seqId';) \n";}
  
  $select_handle->finish;
	$dbh->commit();
  $dbh->disconnect();
  
  return $uniprotId
}


=head2 createEnzymeTable

 Title   : createEnzymeTable
 Usage   : createEnzymeTable( -db => 'databaseName.db',
                                -uniGoPath => '/folder/GODBLink'
			       );

 Function: 		Download the file with the Enzyme Mapping then parses it in such a way that for each EnzymeID are
              taken all the UniprotIDs. If the latter is present then a query INSERT..ON DUPLICATED is performed.
              In this way the Mappings table is complete and an Enzyme table is also created with id and description.
              

 Returns : nothing

=cut
sub createEnzymeTable {
  my $fileToDlPath  = shift;	   
  my ($mappingTableHash) = shift;
  #my $go_enz_create = shift;
  my $uniTable = shift;
  my $uniprotIdHash = shift;#Hash with uniprot Ids mapped with DB ids
  
	my $createEnzTable = "YES";
  my $enzymeTable = 'enzymeTable'; 
  my $enzymeTablePath =  $configHash->{'tempTableFolder'}."/".$enzymeTable;
  
  #This extracts the file name from the web link
  my $fileName = extract_name ($fileToDlPath,"0");
  
  #Download of the file with Enzymes
  if (not (-e $dbDataFolder."/".$fileName) ){
    print "\n Downloading File ".$fileToDlPath."...";
    download_file($fileToDlPath,$dbDataFolder);
    #print "...DONE!\n";
  }else{   print "File ".$fileToDlPath." already downloaded...\n"; }
  
  #Here we create a log file where to write informations on 
  #what are the genes that haven't been found in UNIPROT for the GO terms
  print "If Annocript can't find some Uniprot ID will write in ".$dbDataFolder."/".$configHash->{'uniprotGenesNotPresent'}."\n";
  open (UGLOG, ">>".$dbDataFolder."/".$configHash->{'uniprotGenesNotPresent'}) 
  or print "ERROR: Cannot open log file for missing genes in UNIPROT. Maybe you need more permissions.\n";
     
  #Extract enzyme ids from the file with all the uniprot ids   
  #Here we open the enzyme file from Enzyme Commission and extract a file with ids, descriptions and corresponding Uniprot IDs   
  my $uniprotEnzymeFilePath = extractEnzimeIdsFile($dbDataFolder."/".$fileName);
 
  open(UNIENZYMETABLE,$uniprotEnzymeFilePath);
  open(ENZTABLE, ">$enzymeTablePath");
  my $enzymeIndex = 1;
  my @elements = ();
  my $inserted = 0;
      
  while (my $row = <UNIENZYMETABLE>) {
	  @elements = split("\t",$row); 
    $elements[1] =~ s/[^A-Za-z0-9\-\(\)\s\_\/\.\,]/ /g; #Clean up the description ARE YOU SURE YOU WANT TO CLEAN???
    
    #If the enzymes is with some uniprot ids
    if ( scalar(@elements) > 2){
      $elements[2] =~ s/\s//g;#Erase spaces that are in the string with all the uniprot ids
      my @unipIds = split(",",$elements[2]);
      
      #print_array(\@unipIds);    
      if ($createEnzTable eq "YES"){
       print ENZTABLE join("\t",$enzymeIndex,$elements[0],$elements[1]."\n");
      }
      
      #print "The ids are: ".scalar(@unipIds)."\n";
      foreach my $unipId (@unipIds){
        #Here we search for the uniprot id inside the mapping table (uniprotGO temporarily)
        #We have to add a field at the table and in this field we put the unipIndex
        if ( $uniTable eq 'uniref' ){
               $unipId = getUnirefType($configHash->{'unirefDBLink'}).'_'.$unipId;
        }
         #Usage of a choice by Annocript to use or not an hash that takes a lot of RAM
         #or to do queries to the database     
         my $unipIndex = '';     
         if ( $configHash->{'useHugeRAM'} eq 'NO'){
           $unipIndex = fetchUniprotId($unipId,$uniTable);  
          }elsif ((defined $uniprotIdHash->{$unipId}) or (exists $uniprotIdHash->{$unipId})){
              $unipIndex =  $uniprotIdHash->{$unipId};
          } 

         if ( $unipIndex ne ''){
            #my $unipIndex =  $uniprotIdHash->{$unipId};
            #print "L'id di uniprot per $unipId è: $unipIndex\n";
             $inserted++; 
            if ( (defined $$mappingTableHash->{$unipIndex}) and ( exists $$mappingTableHash->{$unipIndex}->{'e'}) ){
              #Extract all the actual enzyme ids and control if the current is present
              #This control is pratically useless because we are flowing the file of enzymes and they are sorted
              my @terms = split(",",$$mappingTableHash->{$unipIndex}->{'e'});
              pop(@terms);#This is to erase the last element that is simply a space but the split takes it
              if ( scalar(grep { /\b$enzymeIndex\b/} @terms) == 0){
                $$mappingTableHash->{$unipIndex}->{'e'}.=$enzymeIndex.",";
                #print "Adding $enzymeIndex..\n";
              }
                  
            }else{
              $$mappingTableHash->{$unipIndex}->{'e'}.=$enzymeIndex.",";
              #Don't do any CHOP. I'll do later when I write mapping on the tables
              #print "Adding $enzymeIndex..\n";
            }  

        }else {	print UGLOG "WARNING (in createEnzymeTable): there is no $unipId sequence in the uniprotId hash\n";}
       }
    }
    
    $enzymeIndex++;       
  }

  close(UNIENZYMETABLE);
  close(ENZTABLE);
  close(UGLOG);
  
  table_2_db($enzymeTable);
  
  #if ($configHash->{'keepTempFiles'} eq "YES"){
    ##And finally delete the table used that is no more useful
    #delete_file($enzymeTablePath);
    #delete_file($uniprotEnzymeTablePath);
  #}
  
  #delete files downloaded and uncompressed
  if ($configHash->{'eraseFiles'} eq "YES"){
    delete_file($dbDataFolder."/".$fileName);
  }
}

=head2 sortAndWrite

 Title   : sortAndWrite
 Usage   : sortAndWrite( -tableName => 'tableName',
                                -fileName => '/folder/fileForGOUNIPROT'
			       );

 Function: sort the hash table and write it on a file

 Returns : nothing

=cut

sub sortAndWrite{
  my $tableName=shift;
  
  my $path_to_table =  $tempTableFolder."/".$tableName;#Creating the path to the table
  my @sorted_desc = sort {$info->{$tableName}->{$a} <=> $info->{$tableName}->{$b}} keys %{$info->{$tableName}};
  
  open(DESC, ">".$path_to_table);
  foreach my $key (@sorted_desc) {
    #print "key: ".$key. " table: ".$tableName."path_to_table:".$path_to_table."\n";#DEBUGCODE
    print DESC join("\t",$info->{$tableName}->{$key},$key)."\n";
  }
  close(DESC);  
}


=head2 addToHash

 Title   : addToHash
 Usage   : addToHash( -tableName => 'name of the table that have to be written,
                                -tableElement => name ot the element to add
                                -id => id of the element of the table to be incremented
			       );

 Function: adds an element to an hash table. if it still exists its ID is taken

 Returns : nothing

=cut
sub addToHash{
  my $tableName=shift;
  my $tableElement=shift;
  my ($id) = shift;
  
  if( defined($tableElement)){
    #print "$tableName: ".$tableName." ele:  ".$tableElement." id: ".($$id)."\n";#DEBUGCODE
    unless(exists $info->{$tableName}->{$tableElement}){
      $info->{$tableName}->{$tableElement} = ($$id);
      $$id++;
    }
  } 
}



=head2 createDBAndTables

 Title   : createDBAndTables
 Usage   : createDBAndTables( -database => 'name of the database,
                               );

 Function:  creates the database with all the tables

 Returns : nothing

=cut
sub createDBAndTables{
  # CONFIG VARIABLE
	my $undef = undef;
	my $sqlFile= shift;
 
	$sqlFile = $configHash->{'DBCreationFolder'}."/".$sqlFile;
  
  if ( -e $sqlFile){
    #First the database with that name will be erased
    destroy_db($database);

    #print "Creating db ".$database." with user:".$mySqlUser." and pwd: ".$mySqlPass."\n";
    # DATA SOURCE NAME
    my $dsn = "dbi:$platform:information_schema:$host:$port";
    my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } )
              or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
        
    #Then it will be created
    my $rc = $dbh->func("createdb", $database, $host, $mySqlUser, $mySqlPass, 'admin');
     
    if ($dbh->err()) { die "$DBI::errstr\n"; }
    $dbh->disconnect();
    print "Uploading ".$sqlFile." into ".$database."\n"; 
    
    (system qq(mysql -u$mySqlUser -p$mySqlPass $database < $sqlFile) )== 0
    or die "ERROR [$?]: an error occurred while importing $sqlFile in $database. \n"; 
     
    print "DB - proteinDB - and tables created!\n";
    }else {die "ERROR [$?]: the file $sqlFile with the SQL commands to create tables is not present. Annocript will exit: ?\n";}
  
    #Add information about the database in the info table
    my $dbName = $database;
    my $dbRelDate = $configHash->{'dbVersion'};
    my $annocriptVer = $configHash->{'AnnocriptVer'};
    my $dbTablesVersion = $configHash->{'dbTablesVersion'};
    
    query_to_db("INSERT INTO infoTable (dbName, dbRelDate, annocriptVer,dbTablesVersion) ".
                "VALUES ('".$dbName."','".$dbRelDate."','".$annocriptVer."','".$dbTablesVersion."')");
}

=head2 destroy_db

 Title   : destroy_db
 Usage   : destroy_db( -database => 'name of the database,
                               );

 Function:  destroy a database with all the tables

 Returns : nothing

=cut
sub destroy_db{
    my $database = shift;
    
    # DATA SOURCE NAME
    my $dsn = "dbi:$platform:information_schema:$host:$port";
    my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } )
              or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
    
    #First the database with that name will be erased
    $dbh->do("DROP DATABASE IF EXISTS $database") or die "ERROR: unable to cancel database $database! ?"; 
}

=head2 getKBVersion

 Title   : getKBVersion
 Usage   : getKBVersion( - fileAdd = path of the versioning file
                          -folderName = folder where we download the file
                          -user = user for uniprot.org
                          -pass = password for uniprot.org
                               );

 Function:  gets the version of the database file from a file with this information
 Returns : The output is a string with the version 

=cut
sub getKBVersion{

 my $versionDBLink = shift;
 #print "file name: ".$fileAdd."\n";#DEBUGCODE
 my $folderName = shift;
 my $user = shift;
 my $pass = shift;
 my $databaseDl = shift;
 
 my $version = '';
 #my @list = split("/",$versionDBLink);
 my $fileName= extract_name($versionDBLink,0); #pop(@list);
   
 if (not (-e $folderName."/".$fileName) ){
   print "\n Downloading File: ".$versionDBLink."\n";
   download_file($versionDBLink,$folderName, $user,$pass);
 }else{
   print "File ".$folderName."/".$fileName." already downloaded...\n";
 }
 my $line;
  #Open file and read the version
  open(descrFile, $folderName."/".$fileName) or die "ERROR [$?]: Cannot open versioning file...Exiting: ?\n";;
  if ($databaseDl eq "uniprotkb"){ 
   $line = <descrFile>;
   #print "firstLine: ".$firstLine."\n";#DEBUGCODE
   $line =~ m/(^.*)\sKnowledgebase\sRelease\s(.*)\sconsists/;
   #$line =~ m/(^.*)\sRelease\s(.*)\sconsists/;
   $version = $2;
  }
 
  if ($databaseDl eq "uniref"){
    while ($line = <descrFile>){
      #print $line."\n";#DEBUGCODE
      if ($line =~ /Release: (.*),/){
        $version = $1;
      }
    }
  }
  
 if ($version eq ''){
      die "ERROR [$?] Reading the file  of version. Not in the format established in 2013. Please contact the creators of Annocript or ".
      " modify the function getKBVersion in db_creator.pm module. Exiting...: ?\n";
  } 
 close(descrFile);
 
 return $version;
}

1;
