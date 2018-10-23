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
    
 
package LIB::programs_management;

## programs_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of an input configuration file whose variables
#will be put inside an hash

#BEGIN { print (("  " x $main::x++) . "Beginning progmanagement compile\n") }
BEGIN
{
	require Exporter;
	use vars qw(@ISA @EXPORT);
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(  checkConfigVariables configFile2Hash log_and_exit print_and_log 
						initialize_folders checkVariable exec_command try_exec_command);
}
             
use strict;
use warnings;
use Data::Dumper;#To print the hashes
use Cwd;#To change work directory
use Parallel::ForkManager;#For parallel executions
use File::Copy;#To manage files
use LWP::Simple;#For simple URL check (head)

use LIB::files_management qw( );


=head2 exec_command

 Title  : exec_command
 Usage  : exec_command( - command => 'the command to execute',
                      );

 Function: Executes the command given in input by using a sys call

 Returns : nothing

=cut
sub exec_command {
  my $command = shift;
  #print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  ( system ($command)  ) == 0
  or die "ERROR [$?]: an error occurred while running $command: $?";
  #print "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  #print "CORRECTLY DONE command:\n\t$command\n" if $? == 0;
}


=head2 try_exec_command

 Title   : try_exec_command
 Usage   : try_exec_command( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.

 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_command{
    my $command = shift;

    my $maxTimes = 5;
    my $success = -1;
    my $timesCount = 0;

    while ($success == -1 and $timesCount < $maxTimes){
        if ( (system $command) == 0) {
          $success = 1;
        }
        else{
         if ($? == -1) {
              print "failed to execute: $!\n";
          }
          elsif ($? & 127) {
              printf "child died with signal %d, %s coredump\n",
                  ($? & 127),  ($? & 128) ? 'with' : 'without';
          }
          else {
              printf "child exited with value %d\n", $? >> 8;
              
          }
         $timesCount++;
        }
    }
    return $success;
}


=head2 copy_config_to_sess_fold

 Title   : copy_config_to_sess_fold
 Usage   : copy_config_to_sess_fold()
 Function: copies a file into a folder but if it exists needs a specific parameter to force 
			the copy
            
 
 Returns : nothing

=cut
sub copy_file_to_fold{
	my $config_file = shift;
	my $sessionFolder = shift;
	my $force_run  = shift;
	
  print "Copying the given config file in $sessionFolder\n" ;

  #We check if the config file is already present in the given folder. If this happens, the user can decide the overwriting
  if ( -e ($sessionFolder.'/'.$config_file) ){
     exit "WARNING: The file $config_file already exists in $sessionFolder. Please use the parameter --force_run...\n";
  }else{
    (system("cp $ARGV[0] $sessionFolder") ) == 0 
        or annoDie("ERROR: Unable to copy $ARGV[0] in $sessionFolder. Please check permissions...\n");
  }
}



=head2 initialize_folders

 Title   : initialize_folders
 Usage   : initialize_folders(  );

 Function:  this subroutine initializes the basic folders. If the folders.txt file is not present
            means that the program has not been installed in the folder used.
 Returns : nothing

=cut
sub initialize_folders{
  my $foldersFile = shift;

	#print "Opening $foldersFile..\n";
  #If the file with folders references is installed the program can start
  open (FOLD, $foldersFile) or die " ERROR: opening $foldersFile. Please install the program before.\n";
  my $line = <FOLD>;

	#print "Extracting folders path..\n";
  #Extracts folders paths
  my @folders = split(" ", $line);

  #The working folder is the first path
  my $workingFolder = $folders[0];

  #Cancel the final slash..if it is there it is removed
  if($workingFolder =~ /\/$/){
    chop($workingFolder);
  }

  my $programFolder = $folders[1];
  #Cancel the final slash..if it is there it is removed
  if($programFolder =~ /\/$/){
    chop($programFolder);
  }
  close(FOLD);

  return $workingFolder,$programFolder;
}

=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configFile -> file with the user configuration
                              - variablesFile -> the path to a file with all variables written
                              - lineToCheck -> line in the variables file to be used
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.

 Returns : nothing

=cut
sub checkConfigVariables {
  my $configFile = shift;
  my $variablesFile = shift;
  my $lineToCheck = shift;

  my $hashToCheck;

  if (! open(VARF,"<$variablesFile")){ die "ERROR: Failure opening '$variablesFile'. Your program version is corrupted - $!";}
  if (! open(CFILE,"<$configFile")){ die "ERROR: Cannot find '$configFile' - Your program version is corrupted - $!";}

  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CFILE>){
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashToCheck->{$1} = "OK";
    }
  }
  close(CFILE);

  my @confVars = ();

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
      if ($lines == $lineToCheck){
        #print "program Variable: $var - - value: ".$hashToCheck->{$var}."\n";#DEBUGCODE
        #put the variable inside an array for the configuration file
        push (@confVars, $var);
        if( !(defined($hashToCheck->{$var})) ){

          die "ERROR: in $configFile variable $var is missing. Please check the file. Closing...\n ";
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      }
    }
  }

  #print_array(\@allVars);
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashToCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @confVars )){
       die "ERROR: Variable $key is in the config files and not in $variablesFile file. This is completely wrong. Program will not work...\n ";
     }
  }
  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}


=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                             - configHash = the pointer to the hash to be filled
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
			allows to use the empty space to address where the parameter must not be loeaded
 Returns : nothing

=cut
sub configFile2Hash{
  my $configFilePath=shift;
  my ($configHash) = shift;

  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or die "ERROR: The file $configFilePath doesn't exists. The program will exit..\n";
  while (my $line = <configFile>){
		#This IF is useful if you want to put a description above in the text file. Delimit it with a set of hashtags
    if ($line =~ /#########/){
      $start = 1;
    }
   # if( ($line =~ /(\S+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /#/) ){
   if( ($line =~ /(\w+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /^#/) ){
      if ( $2 ne ''){
			  $$configHash->{$1} = $2;
				#annoPrint ("$1 = $2\n") ;#DEBUGCODE
			}
    #}elsif ( $line =~ /(\S+)\s*=/ ){
    }elsif ( $line =~ /(\w+)\s*=/ ){
			delete($$configHash->{$1});	
		}
  }
  close(configFile);
  #print Dumper\$configHash; #DEBUGCODE
}

=head2 checkVariable

 Title   : checkVariable
 Usage   : checkVariable(  - var -> value of the variable to check
                           - name -> name of the variable
                           - sentence -> something to write to the user in case of error)

 Function:  this subroutine checks if a variable has the YES or NO value.
						It takes in input a "yes and a "no sentence" which is used whether 
						the value corresponds

 Returns : nothing

=cut
sub checkVariable {
  my $var = shift;
  my $name = shift;
  my $yesSentence = shift;

  if ( ($var ne 'YES') and ($var ne 'NO')){
    exit "ERROR: Check variable $name in the config file. $var is a wrong value!\n";
  }elsif ($var eq 'YES'){
    print $yesSentence;
  }
}
###############################LOG FUNCTIONS#########################################
=head2 print_and_log

 Title   : print_and_log
 Usage   : print_and_log( - string -> the sentence that have to be print_and_log 
											- onlyLog -> a number);

 Function: will print (the string in input always in the log file and on the STDOUT
					if onlyLog is used then the print will be only in the log file
 Returns : nothing

=cut
sub print_and_log{
  my $string = shift;    
  my $logFile = shift; 
  my $onlyLog = shift;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  if ( defined $onlyLog){
    print LOG $string;
  }else{
    my $STDOUT = *STDOUT;
    my $LOG = *LOG;
    #Prints on both OUT
    for ($LOG, $STDOUT) { print $_ $string; }
  }
  #Close the log
	close(LOG)
}

=head2 log_and_exit

 Title   : log_and_exit
 Usage   : log_and_exit( - string -> the sentence that have to be print in the log

 Function: will print the string and die
 
 Returns : nothing

=cut
sub log_and_exit{
  my $string = shift; 
  my $logFile = shift;     
  
  my $STDOUT = *STDOUT;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  my $LOG = *LOG;
  #Print on both the log and the STDOUT the error
  for ($LOG, $STDOUT) { print $_ $string; }
  #if ( -d $logFolder){
    ##Moving the log file to the log folder before to exit the program
    #move($workingFolder.'/'.$logFile, $logFolder."/".$logFile) or die "Cannot move $workingFolder/$logFile to $logFolder/$logFile\n";
  #}
  #Close the log
	close(LOG);
  exit;
}

1;
   
