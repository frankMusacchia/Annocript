#!/usr/bin/perl

use strict;
use warnings;
use DBI; 
use Switch;#To use the SWITCH construct


  print "1. Check parallel RPS state \n2. Sequences Number of you fasta file \n3. Check BLASTX status \nx. Exit\n";
  my $session = get_input("Your choose: ","^[1234]\$");
    print $session."\n";#Print on the log file
	switch ($session) {
    case ["1"]	{ check_rps_out_state();}
		case ["2"]  { sequences_number();}
    case ["3"]	{ check_blastx_status();}
    case ["x"]  { print "Annocript was interrupted by user... bye!\n";exit;}
    else { die "ERROR: No valid input...\n";}	
  }


#This subroutine says how much sequences contains a fasta file
sub sequences_number {
    my $jobs = 'jobs'; 

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
 my $jobs = 'jobs'; 

 if( (scalar <$jobs/*>) ){
  print "Please insert the name of the folder you want to inspect BLASTx status\n";
   
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

  my $typeIn = get_input( "Please insert the name of the type of BLASTx you want to check [t (Trembl), s (Swiss-Prot), u (Uniref)]\n","^[tsu]\$");
  my $type = '';
   
  switch ($typeIn) {
    case ["t"]	{ $type = "TR";}
		case ["s"]  { $type = "SP";}
    case ["u"]  { $type = "Uniref";}
  }
  
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

# This script will sum the number of sequences already examined by rpsblast in the output folder
# and give the current state
sub check_rps_out_state {
  
 my $fileName = 'rpstblastnOut'; 
 my $jobs = 'jobs'; 

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


#=head2 checkOutFile

 #Title   : checkOutFile
 #Usage   : checkOutFile( -database => 'name of the database,
                               #);

 #Function:  check if an output exists and ask to user if erase or not. Otherwise it does nothing.

 #Returns : nothing

#=cut
#sub checkOutFile{
  #my $outName = shift;
	#if (defined $outName){
  #$outName = $outFolder."/".$outName;
    #if (-e $outName ){
      #my $eraseOut = getInput("\nExistent out file $outName will be erased and all your data lost. Erase it?(y or n)",'^[Y|y|n|N|x|X]$');
      
      #switch ($eraseOut) {
        ##The file exists and the user wants to cancel it
        #case ["y","Y"]	{ my $sure= getInput("Are you sure?(y or n)",'^[YynN]$');
          #switch ($sure) { 
            ##The file exists and the user is sure to cancel it
            #case ["y","Y"] {annoPrint ("File $outName will be overwritten during the execution.\n"); }
            ##The file exists and the user wanted to cancel it (or was wrong in pression of keyboard) and doesn't want anymore
            #case ["n","N"] { print	"Ok, Annocript will use this file...\n";}
            #else { annoDie( "Ok, choose another name in the config file and don't let Annocript waste its time...\n");}					
          #}
        #}
        ##The file exists and the user doesn't wants to cancel it			
        #case ["n","N"] { print	"Ok, Annocript will use this file...\n";}
        #else { annoDie( "Ok, choose another name in the config file and don't let Annocript waste its time...\n");}																
      #}
  #}
  #}
#}



#=head2 findFileByString

 #Title   : findFileByString
 #Usage   : findFileByString( - dir = folder where to search
                          #- whatToSearch  = is the string that indicate what to search
                               #);

 #Function: searches inside a folder if a database of UNIPROT exists by checking if some
            #files with specific extensions exist (.phr,.pin, .psq) and are not empty
			
 #Returns : nothing

#=cut
#sub findFileByString {
	#my $dir = shift;
  #my $whatToSearch = shift;
  
  ##annoPrint ($whatToSearch."\n";
  ##Check if the dir exists
  #opendir DIR, $dir or annoDie( "ERROR : cannot open dir $dir. Check Permissions\n");
  #my @files= readdir DIR;
  #my $dbOk = 0;

  #foreach my $file (@files){
    #if ( ($file =~ /$whatToSearch/) ){
      #$file = $dir."/".$file;
     #if ( !(-z $file) ){
       #$dbOk = 1;
      #}else{$dbOk = 0;last;}
    #}
  #}
  #return $dbOk;
#}

#=head2 table2sql

 #Title   : table2sql
 #Usage   : table2sql( -tableName => 'tableName',
                                #-fileName => '/folder/fileForGOUNIPROT'
			       #);

 #Function: sort the hash table and write it on a file

 #Returns : nothing

#=cut
#sub table2sql{
	#my $table = shift;
	#print "\n Now uploading in db...";
   
  #(system qq(mysqlimport -u$mySqlUser -p$mySqlPass $database $table --local --fields-terminated-by="\t" --lines-terminated-by="\n") )== 0
   #or die "ERROR [$?]: an error occurred while importing $table in $database.\n"; 
    
	#print "....DONE!\n $table Upload completed.\n";
#}
