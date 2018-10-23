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
    
 
package LIB::annocript_management;

## annocript_management.pm
#Author: Francesco Musacchia  (2018)
#Permits the management of an input configuration file whose variables
#will be put inside an hash

#BEGIN { print (("  " x $main::x++) . "Beginning progmanagement compile\n") }
BEGIN
{
	require Exporter;
	use vars qw(@ISA @EXPORT);
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(  check_annocript_integrity valid_analysis check_input_fasta_file
									check_uniref_presence check_db_tables_version check_programs_presence
									 indexing_files_present my_get_proteoms_from_UniProt
									valid_analysis detect_gff_db indexed_db_present getUnirefType
									getSILVADBVersion uniprot_db_selection getDBVersion clean_cd_name
									nice_printing );
}
             
use strict;
use warnings;
use DBI;#DB library
use File::Copy;#To manage file
use Net::FTP;

use LIB::programs_management qw(log_and_exit print_and_log checkVariable );
		
#Using a library to manage files
use LIB::files_management qw( check_presence extract_name);

#Using a library to manage fasta files
use LIB::fasta_management qw( checkFastaFormat detect_fasta_type count_sequences);

#Using a library to manage html and ftp protocols
use LIB::html_management qw(download_file);
								
#Using a library for standard utilities								
use LIB::std_lib qw(correct_type);


#We need these variables to check if the files of the indexed database are present 
#These are for SwissProt, Trembl, UniRef                           
my @blastIndexedExts = qw(.psq .pin .phr);  

#These are for ncRNA.txt
my @ncRNAExts = qw(.nsq .nin .nhr);  



=head2 getSILVADBVersion

 Title   : getSILVADBVersion
 Usage   : getSILVADBVersion( - fileAdd = path of the versioning file
                          -folderName = folder where we download the file
                          -user = user for uniprot.org
                          -pass = password for uniprot.org
                               );

 Function:  gets the version of the SILVA database from a file with this information
 Returns : The output is a string with the version 

=cut
sub getSILVADBVersion{

 my $versionDBLink = shift;  # Link to FTP site of SILVA
    my $version = '';
    my $ftp = Net::FTP->new($versionDBLink) or die "Could not connect to $versionDBLink: $!";
    $ftp->login('anonymous', 'anonymous') or die "Could not login to $versionDBLink: $!";
    $ftp->cwd('current');  # We follow the "current" symbolic link...
    my $current_directory = $ftp->pwd();  # And retrieve the actual name, which looks like "/release_132"
    $current_directory =~ m/release_([0-9]+)/;  # Extract only version number
    $version = $1;
    $ftp->quit();  # Be nice and close the connection
    if ($version eq ''){
        die "ERROR [$?] Reading the file of SILVA version. Not in the format established in 2018. Please contact the creators of Annocript or ".
        "modify the function getSILVADBVersion in db_creator.pm module. Exiting...: ?\n";
    }
    return $version;
}


=head2 getSILVADBVersion

 Title   : getSILVADBVersion
 Usage   : getSILVADBVersion( - fileAdd = path of the versioning file
                          -folderName = folder where we download the file
                          -user = user for uniprot.org
                          -pass = password for uniprot.org
                               );

 Function:  gets the version of the SILVA database from a file with this information
 Returns : The output is a string with the version 

=cut
sub getSILVADBVersionOLD{

 my $versionDBLink = shift;
 my $folderName = shift;
 
 my $version = '';
 my $fileName= extract_name($versionDBLink,0);
   
 if (not (-e $folderName."/".$fileName) ){
   #print "Downloading File: ".$versionDBLink."\n";
   download_file($versionDBLink,$folderName);
 }else{
   print "File ".$folderName."/".$fileName." already downloaded...\n";
 }
 my $line;
  #Open file and read the version
  open(descrFile, $folderName."/".$fileName) or die "ERROR [$?]: Cannot open SILVA version info file...Exiting: ?\n";;
   $line = <descrFile>;
   chop($line);
   #print "firstLine: ".$firstLine."\n";#DEBUGCODE
   my @pieces = split(" ",$line);
   $version = $pieces[scalar(@pieces)-1];
   #$line =~ m/SILVA\sRelease\s(.*)/;
   #$version = $1;
  
 if ($version eq ''){
      die "ERROR [$?] Reading the file of SILVA version. Not in the format established in 2015. Please contact the creators of Annocript or ".
      " modify the function getSILVADBVersion in db_creator.pm module. Exiting...: ?\n";
  } 
 close(descrFile);
 
 return $version;
}

=head2 nice_printing

 Title   : nice_printing
 Usage   : nice_printing( - sentence -> the sentence that have to be annoPrint ();

 Function: this subroutine prints in a box of asterisks the sentence in input
 
 Returns : nothing

=cut
sub nice_printing{
    my $sentence = shift;
    my $log_file = shift;
    
    my $length = length($sentence);
    my $i;
    for ($i=0; $i<($length);$i++){
      print_and_log ("#",$log_file);
    }
    print_and_log ("####\n# $sentence #\n####",$log_file);
    for ( $i=0; $i<$length; $i++){
		print_and_log ("#",$log_file);
    }
    print_and_log ("\n\n",$log_file);
}



=head2 clean_cd_name

 Title   : clean_cd_name
 Usage   : clean_cd_name(  - local folder = the folder where actually is the GFF3_AND_OUTPUT script
                               );

 Function:  gets the hash table with all the path e names in input from the config file in input. Check them and do a lot of other controls
			on input and output files. How they are written, if they are present and so on..
 Returns : nothing

=cut
sub clean_cd_name{
  my $cdName = shift;
   
  #Sometimes in uniprot they put the name of the zip file with LE (Little Endian) and BE (Big Endian)
  if ($cdName =~ /_LE|_BE/){
    $cdName =~ s/_LE|_BE//g;
  }
  return $cdName;
}


=head2 getDBVersion

 Title   : getDBVersion
 Usage   : getDBVersion( - fileAdd = path of the versioning file
                          -folderName = folder where we download the file
                          -user = user for uniprot.org
                          -pass = password for uniprot.org
                               );

 Function:  gets the version of the database file from a file with this information
 Returns : The output is a string with the version 

=cut
sub getDBVersion{

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
   #print "Downloading File: ".$versionDBLink."\n";
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
      die "ERROR [$?] Reading the file  of version. Not in the format expected. Please contact the creators of Annocript or ".
      " modify the function getDBVersion in db_creator.pm module. Exiting...: ?\n";
  } 
 close(descrFile);
 
 return $version;
}


	


=head2 indexed_db_present

 Title   : indexed_db_present
 Usage   : indexed_db_present( -dbToSearch => name of the db that is searching,
                                 -extensions  => an array with the extensions to search,
			       );

 Function: this function checks if a database for the BLAST is consistent

 Returns : returns 1 if the db is ok, otherwise 0

=cut
sub indexed_db_present {
   my $dbToSearch = shift;
   my $chosenExts = shift;#Extensions to search
   my $dbDataFolder = shift;#Folder of dbs


	#First check if Diamond is being used, otherwise check the blast indexing later
	if ($chosenExts eq 'diamondIndexedExts') {  
		# Checking whether DIAMOND .dmnd index exists
    if (-f $dbDataFolder."/".$dbToSearch.".dmnd") { return 1 } else { return 0 }
  }
  # If we're still here, then DIAMOND use it not specified (or BLASTN databases are used which DIAMOND won't support),
  # continue with BLAST check...
    
   #Array with extensions to search ,\@blastIndexedExts,\@ncRNAExts
   my @extensions = ();
      #Needed variables
   my $out;#extension of the file where are indicated the created files
   my $descrWord;#This is a string which specifies when in the previous file
									#starts the list
   my $countFiles = 0;  
   my $shouldBe = 0; 
   my $retVal = 0;
   
   $out  = '.info';
   
   #Get the array with the extensions
   if ($chosenExts eq "blastIndexedExts"){
	    @extensions = @blastIndexedExts;
	    
	    
	}elsif ($chosenExts eq "ncRNAExts"){
		@extensions = @ncRNAExts;
		#$out = ".nal";
		#$descrWord = 'DBLIST';
	}
	$descrWord = 'Volumes:';
  
   #my $dbName = extract_name($dbToSearch,0);
   my $dbName = $dbDataFolder."/".$dbToSearch;
   #go in the dbData folder 
  
   #for the corresponding database we are searching for, open the file with .info extension
   if ( (-e $dbName.$out) and (-s $dbName.$out) ){
     open(VOL_FILE,"<".$dbName.$out) or die "Cannot open file ".$dbName.$out.". Check you permissions.\n"; 
     my $volReached = 0;
     print "Reading the file $dbName$out\n";
     #there you should search for the line where is present "Volumes:"
     while (my $line = <VOL_FILE>){
      if ( $line =~ /$descrWord/ ){
				#If we are searching for uniprot indexing files
				if ($chosenExts eq "blastIndexedExts"){
					while (my $line1 = <VOL_FILE>){
						#then go to the next lines through the end of the file and read all the path it created
						#for each of the lines (and so the paths) you should check that every extension is there
						#to do it we can call the function indexing_files_present        
						chop($line1);#Remove the \n
						$line1 =~ s/\s//g;
						#print "Checking for: $line1...\n";
						#print "db_creator: The extensions are: ";
						#print_array(\@extensions);
						#print "\n";
						
						if (indexing_files_present($line1,\@extensions) == 1){
							 $countFiles +=  scalar(@extensions);
						}
						$shouldBe += scalar(@extensions);
					}
				#If we are searching for the indexed db in the case of nucleotides
				#databases, makeblastdb creates a file with different extension
				}elsif ($chosenExts eq "ncRNAExts"){
					chop($line);#Remove the \n
					my @files = split(" ",$line);
					foreach my $file (@files){
						if ($file ne $descrWord){
							my $filePath = $dbDataFolder."/".$file;
							if (indexing_files_present($filePath,\@extensions) == 1){
								$countFiles +=  scalar(@extensions);
							}
								$shouldBe += scalar(@extensions);
							}
					}
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
		 #If there is no file .out it means there are only three files
		 #so I check those..
		 my $filePath = $dbDataFolder."/".$dbToSearch;
			if (indexing_files_present($filePath,\@extensions) == 1){
							 $countFiles +=  scalar(@extensions);
			}
			$shouldBe += scalar(@extensions);
			
			if ($countFiles == $shouldBe){
       $retVal = 1; 
       print "\n"
      }else{
        print "There should be $shouldBe support files for the indexed database.".
					"  and Annocript found only $countFiles. $dbToSearch must be indexed...\n";
      }  
		 }
   #else{
     #print "Cannot open ".$dbName.$out." The indexing of $dbToSearch will start...\n";
   #}
  return $retVal;
}

=head2 indexed_db_present

 Title   : indexed_db_present
 Usage   : indexed_db_present( -dbToSearch => name of the db that is searching,
                                 -extensions  => an array with the extensions to search,
			       );

 Function: this function checks if a database for the BLAST is consistent

 Returns : returns 1 if the db is ok, otherwise 0

=cut
sub indexed_db_presentOLD {
   my $dbToSearch = shift;
   my $chosenExts = shift;#Extensions to search
   my $dbDataFolder = shift;#Folder of dbs
    
   #Array with extensions to search ,\@blastIndexedExts,\@ncRNAExts
   my @extensions = ();
      #Needed variables
   my $out;#extension of the file where are indicated the created files
   my $descrWord;#This is a string which specifies when in the previous file
									#starts the list
   my $countFiles = 0;  
   my $shouldBe = 0; 
   my $retVal = 0;
   
   $out  = '.info';
   
   #Get the array with the extensions
   if ($chosenExts eq "blastIndexedExts"){
	    @extensions = @blastIndexedExts;
	    
	    
	}elsif ($chosenExts eq "ncRNAExts"){
		@extensions = @ncRNAExts;
		#$out = ".nal";
		#$descrWord = 'DBLIST';
	}
	$descrWord = 'Volumes:';
  
   #my $dbName = extract_name($dbToSearch,0);
   my $dbName = $dbDataFolder."/".$dbToSearch;
   #go in the dbData folder 
  
   #for the corresponding database we are searching for, open the file with .info extension
   if ( (-e $dbName.$out) and (-s $dbName.$out) ){
     open(VOL_FILE,"<".$dbName.$out) or die "Cannot open file ".$dbName.$out.". Check you permissions.\n"; 
     my $volReached = 0;
     #print "Reading the file $dbName$out\n";
     #there you should search for the line where is present "Volumes:"
     while (my $line = <VOL_FILE>){
      if ( $line =~ /$descrWord/ ){
				#If we are searching for uniprot indexing files
				if ($chosenExts eq "blastIndexedExts"){
					while (my $line1 = <VOL_FILE>){
						#then go to the next lines through the end of the file and read all the path it created
						#for each of the lines (and so the paths) you should check that every extension is there
						#to do it we can call the function indexing_files_present        
						chop($line1);#Remove the \n
						$line1 =~ s/\s//g;
						#print "Checking for: $line1...\n";
						#print "db_creator: The extensions are: ";
						#print_array(\@extensions);
						#print "\n";
						
						if (indexing_files_present($line1,\@extensions) == 1){
							 $countFiles +=  scalar(@extensions);
						}
						$shouldBe += scalar(@extensions);
					}
				#If we are searching for the indexed db in the case of nucleotides
				#databases, makeblastdb creates a file with different extension
				}elsif ($chosenExts eq "ncRNAExts"){
					chop($line);#Remove the \n
					my @files = split(" ",$line);
					foreach my $file (@files){
						if ($file ne $descrWord){
							my $filePath = $dbDataFolder."/".$file;
							if (indexing_files_present($filePath,\@extensions) == 1){
								$countFiles +=  scalar(@extensions);
							}
								$shouldBe += scalar(@extensions);
							}
					}
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
		 #If there is no file .out it means there are only three files
		 #so I check those..
		 my $filePath = $dbDataFolder."/".$dbToSearch;
			if (indexing_files_present($filePath,\@extensions) == 1){
							 $countFiles +=  scalar(@extensions);
			}
			$shouldBe += scalar(@extensions);
			
			if ($countFiles == $shouldBe){
       $retVal = 1; 
       print "\n"
      }else{
        print "There should be $shouldBe support files for the indexed database.".
					"  and Annocript found only $countFiles. $dbToSearch must be indexed...\n";
      }  
		 }
   #else{
     #print "Cannot open ".$dbName.$out." The indexing of $dbToSearch will start...\n";
   #}
  return $retVal;
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
  my $configHash = shift;
	my $log_file= shift;
	
	my $ret_val = "NO";
	
  my $dbToSearch = $configHash->{'gffDB'};
	my $platform = $configHash->{'platform'};
	my $host = $configHash->{'host'};
	my $port = $configHash->{'port'};
	my $mySqlUser = $configHash->{'mySqlUser'};
	my $mySqlPass = $configHash->{'mySqlPass'};
  
  # PERL DBI CONNECT AND CREATE TABLE
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } );
  my $databases = $dbh->selectcol_arrayref("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '$dbToSearch'");
  
  my $numKeys = scalar(@$databases);
  
  #The db exists
  if( $numKeys >0 ){
		if ($configHash->{'force_run'}){
			print_and_log("Existent db $dbToSearch. If you chose to obtain the final output, it'll be removed.\n",$log_file);
			$ret_val = "YES";
		}else{
			$ret_val = "NO";
			print_and_log("Existent db $dbToSearch. Choose another name or use --force_run to remove it. Exiting..\n",$log_file);
			exit 1;
		}
    $dbh->disconnect();
  }else{$ret_val = "NO"; }
  
  return $ret_val;
}


=head2 valid_analysis

 Title   : valid_analysis
 Usage   : valid_analysis( $analysis -> the name of the analysis );

 Function:  When user wants to open an existing session Annocript will controll if this does not come
            from a wrong execution: some folders must be inside
 
 Returns : 1 if the analysis is valid, else -1

=cut
sub valid_analysis{
	my $configHash = shift;
  my $analysisFolder = shift;
  
  my $missFolders = "";
  my $retVal = -1;

  if (-d $analysisFolder){		
		#Folders to be alive in the analysis folder if it is consistent
		if ( !(-d $analysisFolder."/".$configHash->{'dataFolder'}) ){
			$missFolders.= " ".$configHash->{'dataFolder'};
		}
		
		if ( !(-d $analysisFolder."/".$configHash->{'logFolder'}) ){
			$missFolders.= " ".$configHash->{'logFolder'}; 
		}
		
		if ( !(-d $analysisFolder."/".$configHash->{'outFolder'}) ){
			$missFolders.= " ".$configHash->{'outFolder'};
		}
		
		if ( !(-d $analysisFolder."/".$configHash->{'statsFolder'}) ){
			$missFolders.= " ".$configHash->{'statsFolder'};
		}
		
		if ( !(-d $analysisFolder."/".$configHash->{'gffFolder'}) ){
			$missFolders.= " ".$configHash->{'gffFolder'};
		}
		
		if ( $missFolders eq ""){
				$retVal = 1;
		}else{
		   print "ERROR: folders $missFolders are missing\n";
		}
 }else{print "Folder $analysisFolder does not exist.\n";}
  return $retVal;
}


=head2 my_get_proteoms_from_UniProt

 Title   : my_get_proteoms_from_UniProt
 Usage   : my_get_proteoms_from_UniProt( - taxon = organism selected
                          -file =  output file
                               );

 Function:  returns a fasta file with peptides for a selected organism in input
					by accessing UniProt directly
 Returns : The output is given as the input fasta file given

=cut
sub my_get_proteoms_from_UniProt {
	my $taxon = shift; # Taxonomy identifier of organism.
	my $file = shift; #Path to the fasta file that will be downloaded
	my $maxTentatives = shift; #Maximum number of tentatives 
	
  open (OUT,">$file");
	my $url = "http://www.uniprot.org/uniprot/?query=organism:$taxon&format=fasta&include=yes";
  my $retVal = 0;
  my $count = 0;
  
  while ( ($retVal == 0) and ($count < $maxTentatives) ){
		print "Tentative $count to query UniProt with: $url...\n"; 
		my $result = get($url);#= $url;
		if (defined $result){
			print OUT $result;
			$retVal = 1;
		}else{
			print "ERROR: accessing $url\n";
			$count++;
		}
	}
  
	close(OUT);
	#print "Sequenze in $file: ";#DEBUGCODE
	#system ( "fgrep -c '>' $file\n");			
	return $retVal;
}


=head2 getUnirefType

 Title   : getUnirefType
 Usage   : getUnirefType( - path = path to the file
                               );

 Function: returns the version of the Uniref database to use as ID: UniRef100, UniRef50 or UniRef90
			
 Returns : a string

=cut
sub getUnirefType {
  my $path = shift;
                  
  my $name = extract_name($path,'1');
  my $type;
  
  if ($name =~ /90/){
    $type = 'UniRef90';
    }
  if ($name =~ /100/){
    $type = 'UniRef100';
    }
  if ($name =~ /50/){
    $type = 'UniRef50';
    }
  
  return $type;
}


=head2 indexing_files_present

 Title   : indexing_files_present
 Usage   : indexing_files_present(  - hash = the reference to the hash to deallocate (\@hash)
                               );

 Function: BLAST indexed databases are composed of a set of files with different extensions.
            If all these files are present we can declare that the database is there
 Returns : nothing

=cut
sub indexing_files_present{
  my $filePath = shift;
  my $extensions = shift;#This needs to be the last parameer since it is an array
  
  my $retVal = 0;
  my $matched = 0;
  
    #print "utilities: The extensions are: ";
    #print_array($extensions);
    #print "\n";  #If for each needed extension the file is there
    foreach my $ext (@$extensions){
      my $fToCheck = $filePath.$ext;
      #print "Check if $fToCheck exists\n";
      if ( (-e $fToCheck) and !(-z $fToCheck) ){
        #print "Wow.. it exists: ".$filePath.$ext."\n";
        $matched++;
      }else{
          print $filePath.$ext." missing or size is zero.\n";
      }
    }
    
  if ( $matched == scalar(@$extensions)){
      $retVal = 1;
  }
  
  return $retVal;
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
	my $configHash = shift;
	my ($blastProgramsInUse) = shift;
	my $log_file = shift;


  #Here we collect the path to the programs using an array initialized in the beginning of the page
  # with all the programs used by Annocript
  foreach my $program (@$blastProgramsInUse){
    my $var= $program.'Path';
    $configHash->{$var} = $configHash->{'blastPath'}.$program;
  }
  
	    #CONTROLS THE PATH TO THE PROGRAMS
    if (!(check_presence($configHash->{'rpstblastnPath'}))){
      log_and_exit($configHash->{'rpstblastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    } 
    
    if (!(check_presence($configHash->{'rpsblastPath'}))){
      log_and_exit($configHash->{'rpsblastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }
    
    if ( (!check_presence($configHash->{'blastxPath'})) ){
      log_and_exit($configHash->{'blastxPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);  
    } 
    
    if ( (!check_presence($configHash->{'blastpPath'})) ){
      log_and_exit($configHash->{'blastxPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);  
    } 
    
    if ( (!check_presence($configHash->{'blastnPath'})) ){
      log_and_exit($configHash->{'blastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }
    
    if ( (!check_presence($configHash->{'tblastnPath'})) ){
      log_and_exit($configHash->{'blastnPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }
    
    if ( (!check_presence($configHash->{'makeblastdbPath'})) ){
      log_and_exit($configHash->{'makeblastdbPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }
    
    if ( (!check_presence($configHash->{'lncRNAPredProgPath'})) ){
      log_and_exit($configHash->{'lncRNAPredProgPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }
    
    if ( (!check_presence($configHash->{'dna2pepPath'}))  ) {
      log_and_exit($configHash->{'dna2pepPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    }

		##Check diamond presence
		#if ( (!check_presence($configHash->{'diamondPath'}))  ) {
      #log_and_exit($configHash->{'dna2pepPath'}." is missing. Please check the links in the configuration file or re-install the program. Exiting...\n",$log_file);
    #}
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
	my $configHash = shift;
	my $log_file = shift;

  my $dbToSearch = $configHash->{'database'};
  my $platform = $configHash->{'platform'};
  my $host = $configHash->{'host'};
  my $port = $configHash->{'port'};
  my $mySqlUser = $configHash->{'mySqlUser'};
  my $mySqlPass = $configHash->{'mySqlPass'};
  my $doDbCreation = $configHash->{'doDbCreation'};
  my $doBuildOutput = $configHash->{'doBuildOutput'};
  my $dbTablesVersion = $configHash->{'dbTablesVersion'};
  
   
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
      check_db_tables_version($dsn2,$mySqlUser,$mySqlPass,$dbTablesVersion) 
				or log_and_exit( "In this version of Annocript (".$configHash->{'AnnocriptVer'}.") the database will have different structure.". 
													"You cannot update a previous version. Please create a completely new database.\n",$log_file);
      if ($configHash->{'force_run'} eq ''){
				log_and_exit( "WARNING: The database $dbToSearch already exists. It can be checked and completed. Use the parameter --force_run. Exiting...\n",$log_file); 
			}else{
				print_and_log( "The database $dbToSearch already exists. It will be checked and completed.\n",$log_file); 
				$configHash->{'dbExistent'} = 'YES';
			}
      $dbh->disconnect();
    }else{$configHash->{'dbExistent'} = 'NO';}
  }
  #Do not create database for proteins but it already exists and we want to use to build output
  elsif(( $doDbCreation eq "NO") and ( $doBuildOutput eq "YES") ){
    #If it doesn't exists than user has to reconfigure Annocript
    if( $numKeys == 0 ){
      print_and_log( "ERROR: Database $dbToSearch you choosed in the config file doesn't exists. Please check it. Annocript will exit...\n",$log_file);
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
  my $dbTablesVersion = shift;
  
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




=head2 check_input_fasta_file

 Title   : check_input_fasta_file
 Usage   : check_input_fasta_file(  - fastaSeqs = the name of the fasta
                                                                                                                                       $
                                                                                                                                       $
                                                                                                                                       $
                                                                                                                                       $
                               );

 Function:   check if the query FASTA file is well written and is in the working folder. Behaviour of the function is different
          if user is calling a saved session. The file is moved soon in the session folder.


 Returns : In the config hash will be returned:
                                               	- $$configHash->{'Present'} => YES if the file is ok, NO otherwise
                                                - $$configHash->{'sequencesAlphabet'} => alphabet type
                                                - $$configHash->{'totSequences'} => total number of sequences in the file

=cut
sub check_input_fasta_file {

#CONTROLS IF THE QUERY EXISTS AND THE NAME IS .fa or .fasta
  my $fastaSeqs = shift;
  my ($configHash) = shift;#An hash with variables to be used in the main program
  my $newSession = shift;#Indicates if the session is loaded
  my $workingFolder = shift;#The folder where the user is working
  my $sessionFolder = shift;#The folder used for the output

  my $name = '';

  #The name has to be without dots inside. Only one dot can separate the extension from the name
  my $dots=0;
  my @line = split (//,$fastaSeqs);
  foreach my $char(@line){
      if ($char eq '.'){
          $dots++;
      }
  }
  if ($dots>1){
      die "Please change the name of your file removing all internal dots (not the one used for the extension). Exiting...\n";
  }
                                                                                                                                        
  #annoPrint ("Fasta name read: ".$fastaSeqs."\n";
  my @name = split(/\./,$fastaSeqs);
  print  "Fasta name read: ".$name[0]."\n";#DEBUGCODE
  #Here is created the name to use for all the files of output by extracting it from the fasta file name
  #It cannot be longer than 50 chars
  if (@name > 1){
    #$$configHash->{'name'} = extract_name($$configHash->{'sessionFolder'},0);
    $name = $name[0];
    #Checks for permitted characters
    if ( (length($name) > 50) # or !($$configHash->{'fastaSeqs'} =~ /(\.fasta?)|(\.fa?)/)
            or ($fastaSeqs !~ /^[A-za-z0-9\_\-]+\.(fa|fasta)$/i)){
            die "$fastaSeqs is not a correct FASTA file name. Allowed characters [A-Z,a-z,0-9,_,-]. Allowed extensions [.fasta,.fa]. Max length: 50 chars\n";
    }else{
      #If session is a new one
      #Fasta file must stay only in the working folder and it will be shifted. If it is not there,then only db creation can be done or $
      if ( $newSession ){
       my $fastaPath = $workingFolder.'/'.$fastaSeqs;
       print "Path to control for fasta: $fastaPath\n";
        #The Sequences file can also stay already in the session folder created by the user in a previous execution
        if (-e $fastaPath) {
          print "Sequences found in your working folder. Checking...\n";
          #After the name all the fasta is checked to see if it respects the standards
			my $warnings = checkFastaFormat($fastaPath,$$configHash->{'maxFastaSeqLen'},$$configHash->{'nuclIUPAC'},$$configHash->{'protIUPAC'});
			if($warnings > 0){
					die "Fasta check failed! The transcriptome you are using have the listed problems ".
					" which can cause stochastic error in BLAST programs execution. ".
					"Please try to correct the problems. Exiting...\n";
		   }
          print "Copying the fasta in session folder...\n";
          #move($fastaPath,$sessionFolder);
          copy($fastaPath,$sessionFolder) or annoDie("Cannot copy $fastaPath in $sessionFolder");
          $$configHash->{'fastaSeqsPresent'} = 'YES';
        }
        else{
          #This else happens when sequences file is not there. Then, if the user chose to execute analysis or to print output
          # the program has to die.
          if ( ($$configHash->{'doExecutePrograms'} eq 'YES') or ($$configHash->{'doBuildOutput'} eq 'YES')
            or ($$configHash->{'extractStatistics'} eq 'YES')){
              die "No sequences file found in $fastaPath. You'll cannot execute analysis...\n";
            }
          print "No sequences file found in  $fastaPath. But you are only creating a DB...\n";
          $$configHash->{'fastaSeqsPresent'} = 'NO';
        }
      }else{#If session is loaded then the fasta must stay only in the session folder. No other fasta can be shifted there. Die otherwi$
        if (-e $sessionFolder.'/'.$fastaSeqs) {
          print "Sequences found in your session folder. Checking it...\n";
          $$configHash->{'fastaSeqsPresent'} = 'YES';
                                #After the name all the fasta is checked to see if it respects the standards
					my $warnings = checkFastaFormat($sessionFolder.'/'.$$configHash->{'fastaSeqs'},$$configHash->{'maxFastaSeqLen'},$$configHash->{'nuclIUPAC'},$$configHash->{'protIUPAC'});
					if($warnings > 0){
									die "Fasta check failed! The transcriptome you are using have the listed problems ".
									" which can cause stochastic error in BLAST programs execution. ".
									"Please try to correct the problems. Exiting...\n";
					}
        }else{
          die "Unable to find ".$$configHash->{'fastaSeqs'}." in $sessionFolder. The program will exit and you can figure out why.".
                    " \nThere are 2 possible reasons: \n\t 1. you written a wrong fasta name in the configuration file.".
                    "\n\t 2. the name in the config file is right but not the folder you chose\n Bye!\n ";
        }
      }

      if ($$configHash->{'fastaSeqsPresent'} eq 'YES'){

        #Here we control what kind of alphabet has the sequence and store the value in a variable
        $$configHash->{'sequencesAlphabet'} = detect_fasta_type($sessionFolder.'/'.$fastaSeqs);
        print "The program detected that your sequences are ";
        if ( ($$configHash->{'sequencesAlphabet'} eq 'dna') or ($$configHash->{'sequencesAlphabet'} eq 'rna') ){
            print "transcripts\n";
        }else{
            print $$configHash->{'sequencesAlphabet'}."\n";
         }
        #Here we count and store the number of sequences in the fasta file
        my $seqsPath = $sessionFolder.'/'.$fastaSeqs;
        my $totSequences = count_sequences($seqsPath); #Extract number of sequences from the query file
        print "Number of sequences: $totSequences";
        $$configHash->{'totSequences'} = $totSequences;
      }
    }
  }else {
    die "The name of file with transcripts needs to have a .fa or .fasta extension. Exiting...\n";
  }
  
  return $name;
}


=head2 check_annocript_integrity

 Title   : check_annocript_integrity
 Usage   : check_annocript_integrity( - configHash = the piece of the config hash for PROGRAMS_EXEC
							 - local folder = the folder where actually is the PROGRAMS_EXEC script
                               );

 Function:  checks variables mostly inside the annocript configuration file. 
 Returns : nothing

=cut
sub check_annocript_integrity {
	my $configHash  =shift;
	my $log_file = shift;
 
  #######CHECKING integrity of Annocript
  #print_and_log("Checking folders integrity...\n",$log_file);#DEBUGCODE
  if ( ! ( (-e $configHash->{'lib_folder'}."/".$configHash->{'R_seqs_stats_script'})
	   and (-e $configHash->{'lib_folder'}."/".$configHash->{'R_piechart_script'})
	   and (-e $configHash->{'lib_folder'}."/".$configHash->{'R_barplot_script'})
        and (-d $configHash->{'lib_folder'}) 
              and (-d $configHash->{'usefulFolder'}) 
                and (-d $configHash->{'guideFolder'}) 
                  and (-d $configHash->{'configFolder'}) ) ){      
    log_and_exit( "ERROR: This Annocript version is corrupted. Please download and re-install Annocript. Exiting..\n",$log_file);
  }    
  
  #Checking folders and files names
  #print_and_log("Checking folders and files names...\n",$log_file);#DEBUGCODE
  #Names of the output of the programs
  if ( $configHash->{'rpstblastnOut'} !~ /^[A-za-z0-9\_\-]*$/i  ){
		log_and_exit( "The name ".$configHash->{'rpstblastnOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'blastxSPOut'} !~ /^[A-za-z0-9\_\-]*$/i  ){ 
		log_and_exit( "The name ".$configHash->{'blastxSPOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'blastxTROut'} !~ /^[A-za-z0-9\_\-]*$/i  ){
		log_and_exit( "The name ".$configHash->{'blastxTROut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'blastxUnirefOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name".$configHash->{'blastxUnirefOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'blastnOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'blastnOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'lncRNAPredOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'lncRNAPredOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'dna2pepOut'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'dna2pepOut'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'outCount'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'outCount'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  #Names of the stats files
  #print_and_log("Checking Names of the stats files...\n",$log_file);#DEBUGCODE
  if ( $configHash->{'bpStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'bpStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'mfStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'mfStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'ccStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'ccStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'cdStatsFile'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'cdStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
    if ( $configHash->{'pathwaysL1File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'mfStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'pathwaysL2File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'ccStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'pathwaysL3File'} !~ /^[A-za-z0-9\_\-]*$/i ){
		log_and_exit( "The name ".$configHash->{'cdStatsFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  #Filenames with a dot
  # print_and_log("Checking Filenames with a dot...\n",$log_file);#DEBUGCODE
 
  if ( $configHash->{'outFileName'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'outFileName'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'outFiltered'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'outFiltered'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'ORFFastaFileName'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'ORFFastaFileName'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'NCOutFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'NCOutFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'gcContentFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'gcContentFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'closerOSFile1'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'closerOSFile1'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
	if ( $configHash->{'closerOSFile2'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'closerOSFile2'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  if ( $configHash->{'ATGCPercFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'ATGCPercFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  if ( $configHash->{'outHashFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'outHashFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
 
  if ( $configHash->{'parsingFilePath'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'parsingFilePath'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  if ( $configHash->{'filtMappingFile'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'filtMappingFile'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  if ( $configHash->{'R_barplot_script'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'R_barplot_script'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}  

  if ( $configHash->{'R_piechart_script'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'R_piechart_script'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}  
	  
  if ( $configHash->{'R_seqs_stats_script'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'R_seqs_stats_script'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
  }  
	
  if ( $configHash->{'R_log_file'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'R_log_file'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  if ( $configHash->{'gffDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'gffDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	} 
  
  #if ( $configHash->{'ourRiboSeqs'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		#log_and_exit( "The name".$configHash->{'ourRiboSeqs'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file));
	#}  
  
  if ( $configHash->{'ncDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'ncDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}
  
  if ( $configHash->{'gffDB'} !~ /^[A-za-z0-9\_\-\.]*$/i ){
		log_and_exit( "The name ".$configHash->{'gffDB'}." isn't permitted. Choose another in the config file using [A-za-z0-9\_\-]. Program will exit...\n",$log_file);
	}   
  
  #Filenames with HTML extension
   #  print_and_log("Checking Filenames with HTML extension...\n",$log_file);#DEBUGCODE

  if ( ($configHash->{'htmlHome'} !~ /^[A-za-z0-9\_\-\.]*$/i ) and ($configHash->{'htmlHome'} !~ /.html/ )){
		log_and_exit( "The name ".$configHash->{'htmlHome'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n",$log_file);
  }   
  if ( ($configHash->{'htmlAnnotStats'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'htmlAnnotStats'} !~ /.html/ ) ){
		log_and_exit( "The name ".$configHash->{'htmlAnnotStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n",$log_file);
	}   
  if ( ($configHash->{'htmlSeqStats'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'htmlSeqStats'} !~ /.html/ ) ){
		log_and_exit( "The name ".$configHash->{'htmlSeqStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the .html extension. Program will exit...\n",$log_file);
  }
	
 #Logo name check
  if ( ($configHash->{'AnnocriptLogo'} !~ /^[A-za-z0-9\_\-\.]*$/i) and ($configHash->{'AnnocriptLogo'} !~ /.jpeg|.png/ ) ){
		log_and_exit( "What will you use as logo?! ".$configHash->{'htmlSeqStats'}." isn't permitted. Use only [A-za-z0-9\_\-] and the png or jpeg format. Program will exit...\n",$log_file);
  }
		  
  #Checking other variables 
  # print_and_log("Checking other variables ...\n",$log_file);#DEBUGCODE
 
  checkVariable($configHash->{'keepTempFiles'},'keepTempFiles',"Temporary files will be stored\n",$log_file);
 
  if (!correct_type($configHash->{'maxLengthDescs'},"positiveint")){
    log_and_exit( $configHash->{'maxLengthDescs'}." is not a positive integer. Please change it and re-execute Annocript. Exiting...\n",$log_file);
  } 
  
  if (!correct_type($configHash->{'bigMemory'},"positiveint")){
    log_and_exit( $configHash->{'bigMemory'}." is not a positive integer. Please change it and re-execute Annocript. Exiting...\n",$log_file);
  }
}


1;
