#############################utils.pm####################
# This module is to execute few utilities in Annocript  #
#########################################################

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

package USEFUL::utils;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( execute_utils );


use strict;
use warnings;
use DBI; 
use Switch;#To use the SWITCH construct
use Cwd;#To change work directory

#Annocript libraries
use USEFUL::utilities qw( extract_name download_file try_exec_command);
                  
#Global variables
my $configHash;#This hash is to keep variables from Annocript

#Variables present also in annocript.pl (IF YOU CHANGE THERE, CHANGE HERE and VICEVERSA)
my $jobs = 'jobs';#jobs folder
my $mainDataFolder = "data";#Folder inside Annocript with all data needed more than one time
my $annocriptConfig;

#Main subroutine
sub execute_utils {
 $annocriptConfig = shift;
 print "Annocript configuration file: ".$annocriptConfig."\n";	
 configFile2Hash($annocriptConfig);#Reads the configuration file of Annocript and stores the variables
 my $exit = 0;
 while ($exit == 0){
  print "1. Sequences number of you fasta file\n2. Check BLASTX status \n3. Check parallel RPS state ".
  "\n4. Check Portrait state \n5. Search organism\n6. Exit\n";
  my $session = get_input("Your choose: ","^[123456]\$");
	print $session."\n";#Print on the log file
	if ($session eq '1' ) { sequences_number();}
	elsif ($session eq '2' )	{ check_blastx_status();}
	elsif ($session eq '3' )	{ check_rps_out_state();}
	elsif ($session eq '4' )	{ check_portrait_status();}
	elsif ($session eq '5' )	{ search_organism();}
	elsif ($session eq '6' ) { 
		$exit = 1;
		print "Thank you for using Annocript... bye!\n";
		exit;
	}
   else { die "ERROR: No valid input...\n";}	
 }
}

#This subroutine downloads a file with all the organisms of UniProt
sub search_organism {
	my $orgListLink = $configHash->{'uniprotOrgListLink'};
	
	#print "Link: ".$orgListLink;
	#This extracts the file name from the web link
    my $orgListFileName = extract_name ($orgListLink,"0");
    
    #Download of the file with organisms list
    if (not (-e $mainDataFolder."/".$orgListFileName) ){############
      print "\n Downloading File ".$orgListLink."...";
      download_file($orgListLink,$mainDataFolder);
      #print "...DONE!\n";
    }else{   print "File ".$orgListLink." already downloaded...\n"; }
      
    print "A file with all the organisms available in UniProt has been downloaded. Do you want to:\n";
    
    my $exit = 0;
        
   while ($exit == 0){
		my $choose = get_input( "\n1. Read list of organisms\n2. Search Organism\n3. Exit\nYour Choice: ","^[123]\$");
		
		if ($choose == 1){
			(system("more ".$mainDataFolder."/".$orgListFileName)) == 0 
					or die "Can't open file ".$mainDataFolder."/".$orgListFileName."...\n";  
		}elsif($choose == 2){
			print "\nOrganism to search: \n";            
			my $organism = <STDIN>;
			chomp $organism;
			(system("grep 'N=' ".$mainDataFolder."/".$orgListFileName." | grep -i '$organism' " )) == 0 
					or print "Cannot find '$organism' in ".$mainDataFolder."/".$orgListFileName."...\n";  
		}elsif($choose == 3){
			$exit = 1;
		}
	}
}

#This subroutine says how much sequences contains a fasta file
sub sequences_number {

 if( (scalar <$jobs/*>) ){
    print "Please insert the name of the folder you want to inspect number of sequences\n";
     
    (system("ls $jobs/")) == 0 
                or die "Cant open folder $jobs...\n";
    print "\nYour choose: ";            
    my $myFolder = <STDIN>;
    chomp $myFolder;
    $myFolder = $jobs."/".$myFolder;
    
    while( !(-d $myFolder) ){
      (system("ls $jobs/") ) == 0 
      or die "Unable to list directory $jobs!\n";
      print "Choose your folder: ";
      $myFolder = <STDIN>;
      chomp $myFolder;
      $myFolder = $jobs."/".$myFolder;
    }
    
     #print "DEBUGCODE: Entered since the folder exists..\n";
		opendir DIR, $myFolder or die "ERROR [$?]: cannot open dir $myFolder \n";
		my @inFolder = readdir DIR;
		closedir DIR;
    if (scalar( @inFolder)){
      foreach my $file (@inFolder){
        if ( ( ($file =~ 'fasta') or ($file =~ 'fa') )   and (!($file =~ 'portrait'))  and (!($file =~ 'cleanh')) ){
           my $pathFile = $myFolder."/".$file;
            print "Number of sequences in $file: ";
           print `fgrep -c '>' $pathFile`;
           print "\n";
        }
      }
    }
  }
}


# This script will sum the number of sequences already examined by blastx in the output folder
# and gives the current status
sub check_blastx_status {
  
 my $fileName = 'blastx'; 

 if( (scalar <$jobs/*>) ){
  print "Please insert the name of the folder you want to inspect BLASTx status\n";
  my $workDir = getcwd;
  chdir $jobs;#go to jobs directory
  (system("ls -d */")) == 0 
              or die "Cant open folder $jobs...\n";
  print "\nYour choose: ";            
  my $myFolder = <STDIN>;
  chomp $myFolder;
  
  while( !(-d $myFolder) ){
    (system("ls -d */") ) == 0 
    or die "Unable to list directory $jobs!\n";
    print "Choose your folder: ";
    $myFolder = <STDIN>;
    chomp $myFolder;

  }
  chdir $workDir;#back to the working directory
  
  
  my $typeIn = get_input( "Please insert the name of the type of BLASTx you want to check [t (Trembl), s (Swiss-Prot), u (Uniref)]\n","^[tsu]\$");
  my $type = '';
   
  switch ($typeIn) {
    case ["t"]	{ $type = "TR";}
		case ["s"]  { $type = "SP";}
    case ["u"]  { $type = "Uniref";}
  }
  $myFolder = $jobs."/".$myFolder;
  #my $totSequences  = `fgrep -c '>' $myFolder`;
  #Now entering in the rps out folder and checking if it is empty
  my $blastxOutFolder = $myFolder.'/output/blastx';
  if( (scalar <$blastxOutFolder/*>) ){
    my $totDone = 0;
    
    #print "DEBUGCODE: Entered since the folder exists..\n";
		opendir DIR, $blastxOutFolder or die "ERROR [$?]: cannot open dir $blastxOutFolder \n";
		my @blastOuts= readdir DIR;
		closedir DIR;
    
    my $numOutFiles = scalar(@blastOuts)-2;#This is necessary to not count . and .. folders
		if ($numOutFiles >0 ){
      print "\nYou have $numOutFiles output files inside the folder $blastxOutFolder...\n ";
      my $index;
			for ($index = 1; $index <= $numOutFiles; $index++){
				my $blastOut = $blastxOutFolder."/".$fileName.$type."Out".$index;
				 #print "DEBUGCODE: trying to access $blastOut\n";
				#Here we access the blast out only if it corresponds to the specific algorithm
				if ( (-e $blastOut) ){#&& ($blastOut =~ /$fileName/)){
            my $done = `fgrep -c 'Query=' $blastOut`;
            #print "done: $done - ";
            $totDone = $totDone + $done;
        }else{
            print "No file $blastOut has been found..\n ";
        }
      }
    }
   print "Sequences blastxed for $fileName$type=> $totDone\n"; 
  }else {
    die "The folder $blastxOutFolder is empty\n";
    }
  
 }else{
    die "The folder $jobs is empty\n";
  } 
}

# This script will sum the number of sequences already examined by blastx in the output folder
# and gives the current status
sub check_portrait_status {
  
 my $fileName = 'portraitOut'; 

 if( (scalar <$jobs/*>) ){
  print "Please insert the name of the folder you want to inspect Portrait status\n";
  my $workDir = getcwd;
  chdir $jobs;#go to jobs directory
  (system("ls -d */")) == 0 or die "Cant open folder $jobs...\n";
  print "\nYour choose: ";            
  my $myFolder = <STDIN>;
  chomp $myFolder;
  
  while( !(-d $myFolder) ){
    (system("ls -d */") ) == 0 or die "Unable to list directory $jobs!\n";
    print "Choose your folder: ";
    $myFolder = <STDIN>;
    chomp $myFolder;
  }
  chdir $workDir;#back to the working directory

  $myFolder = $jobs."/".$myFolder;
  #Now entering in the rps out folder and checking if it is empty
  my $portraitOut = $myFolder.'/output/'.$fileName;
  my $done;
  if ( (-e $portraitOut) ){
		my $command = "wc -l $portraitOut";
    try_exec_command($command) > 0 or die "Error: unable to execute $command\n";
  }else{
    print "No file $portraitOut has been found..\n ";
  }
 }else{
    die "The folder $jobs is empty\n";
 } 
}

# This script will sum the number of sequences already examined by rpsblast in the output folder
# and give the current state
sub check_rps_out_state {
  
 my $fileName = 'rpstblastnOut'; 


 if( (scalar <$jobs/*>) ){
  print "Please insert the name of the folder you want to inspect RPS state\n";
   
  (system("ls $jobs/")) == 0 
              or die "Cant open folder $jobs...\n";
  print "\nYour choose: ";            
  my $myFolder = <STDIN>;
  chomp $myFolder;
  $myFolder = $jobs."/".$myFolder;
  
  while( !(-d $myFolder) ){
    (system("ls $jobs/") ) == 0 
    or die "Unable to list directory $jobs!\n";
    print "Choose your folder: ";
    $myFolder = <STDIN>;
    chomp $myFolder;
    $myFolder = $jobs."/".$myFolder;
  }
  
  
  #my $totSequences  = `fgrep -c '>' $myFolder`;
  #Now entering in the rps out folder and checking if it is empty
  my $rpsOutFolder = $myFolder.'/output/rpstblastn';
  if( (scalar <$rpsOutFolder/*>) ){
    my $totDone = 0;
    
    #print "DEBUGCODE: Entered since the folder exists..\n";
		opendir DIR, $rpsOutFolder or die "ERROR [$?]: cannot open dir $rpsOutFolder \n";
		my @blastOuts= readdir DIR;
		closedir DIR;
    
    my $numOutFiles = scalar(@blastOuts)-2;#This is necessary to not count . and .. folders
		if ($numOutFiles >0 ){
      print "\nYou have $numOutFiles output files inside the folder $rpsOutFolder...\n ";
      my $index;
			for ($index = 1; $index <= $numOutFiles; $index++){
				my $blastOut = $rpsOutFolder."/".$fileName.$index;
				 #print "DEBUGCODE: trying to access $blastOut\n";
				#Here we access the blast out only if it corresponds to the specific algorithm
				if ( (-e $blastOut) ){#&& ($blastOut =~ /$fileName/)){
            my $done = `fgrep -c 'Query=' $blastOut`;
            #print "done: $done - ";
            $totDone = $totDone + $done;
        }
      }
    }
   print "Sequences rpsblasted => $totDone\n"; 
  }else {
    die "The folder $rpsOutFolder is empty\n";
    }
  
 }else{
    die "The folder $jobs is empty\n";
  } 
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
    print $sentence." ";
    $input = <STDIN>;
    chomp $input;
    print "\n";
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

1;
