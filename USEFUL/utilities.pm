#######################utilities.pm###########################################################
# utilities - A module that contains a series of utils subroutine                              #
################################################################################################

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

package USEFUL::utilities;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( checkLink correct_type db_present delete_file dl_and_extract 
                  download_file extract_file my_extract_file extract_name detect_fasta_type 
                  append_file_2_file clean_cd_name print_array delete_directory 
                  num_processors exec_command count_sequences print_file check_presence
                  my_head try_exec_command deallocate_hash indexing_files_present ram_memory
                  is_folder_empty getUnirefType save_hash join_files_with_cat
                  extract_special_chars create_slices);

#GENERAL UTILITY MODULES
use strict;
use warnings;
use Switch; #To use the SWITCH statement

#FILE MANAGEMENT MODULES
use File::Path qw(make_path remove_tree ) ;
use File::Fetch; #TO use FTP download functions

#HTTP FTP CHECKS
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Response;

##################GENERIC OPERATIONS ON FILES


=head2 delete_file

 Title   : delete_file
 Usage   : delete_file( -filePath => 'a path to a file to delete'
			       );

 Function:  Delete a file in a given location. It can erase also more than one file starting with the same name. This
            can be done by using the usual '*'. When it is found the 'glob' keyword will be used.
 Returns : nothing

=cut
sub delete_file{
  my $filePath = shift;
  
  
    #With this first IF we can delete a set of file with the same name
    if ($filePath =~ /\*/){
      unlink glob $filePath or warn "$filePath  was not deleted for some problems with unlink subroutine or file don't exists anymore.\n";
    }elsif ( -e $filePath ){
        if  ( unlink($filePath) == 1) {
          print "$filePath deleted successfully.\n";
        }else{
          print "$filePath  was not deleted for some problems with unlink subroutine or file don't exists anymore.\n";
        }
    }else{
        print "$filePath has been already erased.\n";
    }
}

=head2 delete_directory

 Title   : delete_directory
 Usage   : delete_directory( -filePath => 'a path to a directory to delete'
			       );

 Function:  Delete a directory with all its content in a given location
 Returns : nothing

=cut
sub delete_directory{
  my $dirPath = shift;
  
  remove_tree($dirPath);
}  


=head2 dl_and_extract

 Title   : dl_and_extract
 Usage   : dl_and_extract( -link = the link to the file on the web
                            - compressed = the file name compressed
                             - uncompressed = the file name uncompressed
                              - user = the username to access some ftp (optional)
                               - pass = the password to access some ftp (optional)
                               );

 Function:  takes in input a link to a compressed file and downloads it if the uncompressed o the compressed
            are not in the folder. The user and password informations need only in some cases. So you can
            give to the function these two last parameters too.

 Returns : nothing


=cut
sub dl_and_extract{
  my $link = shift;
  my $compressed = shift;
  my $uncompressed = shift;
  my $dbDataFolder = shift;
  my $user = shift;
  my $pass = shift;
  #my $originalFileDim = shift;
  
  #Download the file 
 unless (-e $dbDataFolder."/".$uncompressed){#Doesn't exists the  file, proceed
		unless (-e $dbDataFolder."/".$compressed) {#Doesn't exists the compressed, download it
			 print "\n Downloading File: ".$link."\n";
        if (defined ($user) and defined ($pass)){
           download_file($link,$dbDataFolder,$user,$pass);
        }else{
           download_file($link,$dbDataFolder);
        }
        ##Check if the dimensions are  equal to those taken during the link checking
        ##If original is 0 it means that the check has been incorrect
        #if ( $originalFileDim > 0){
          #if (compare_dimensions($originalFileDim,$dbDataFolder."/".$compressed) <= 0){
              #die "Error: ".$dbDataFolder."/".$compressed." has not been well downloaded. Delete it before to start again Annocript.\n";
              ##delete_file($dbDataFolder."/".$compressed);
          #}
        #}
		 }else{ print "File ".$dbDataFolder."/".$compressed." already downloaded...\n"; }
		 print "\n Uncompress File: ".$dbDataFolder."/".$compressed."-->".$dbDataFolder."/".$uncompressed."\n";
    my_extract_file($dbDataFolder."/".$compressed,$dbDataFolder);
 }else{print "File ".$dbDataFolder."/".$uncompressed." already present...\n";}#Exists the file 
}



=head2 download_file

 Title   : download_file
 Usage   : download_file( -fileAdd => is the path of the file
                          -folderName => is the folder where we have to download
                          -user => username froruniprot
                            -pass => password for uniprot
                               );

 Function:  Given in input a path to a file, it download it in the data folder. It assigns as name of the file the last
#element of the array composed by the split of the string path

 Returns : nothing

=cut
sub download_file{
  my $fileAdd = shift;
  #print "file name: ".$fileAdd."\n";
  my $folderName = shift;
  my $user = shift;
  my $pass =shift;
  my $originalFileDim = shift;
  
  
  my $ua = LWP::UserAgent->new;
  my $req = HTTP::Request->new (GET => $fileAdd);
  if (defined($user) and defined ($pass)){
    $req->authorization_basic($user, $pass);
  }

  my $fileToDl = extract_name($fileAdd,"0");
  my $res = $ua->request($req, $folderName."/".$fileToDl);

  if ($res->is_success) {
     print $fileAdd." downloaded successfully!\n";
      
  }
  else {
	  #Here we add a specific control and re-download because sometimes for problems of networks
	  #the file is not downloadable but it becomes after some seconds.
     print "Some problems occurred while downloading $fileAdd. Trying again...\n";
	   my $res2 = $ua->request($req, $folderName."/".$fileToDl);
	   if ($res2->is_success) {
			print $fileAdd." downloaded successfully!\n";
      
      }else{
     
        #This command will force the download with different methods
        if (ftp_fetch_file($fileAdd,$folderName,$user,$pass) == 1){
            print $fileAdd." downloaded successfully!\n";
        }else{
          delete_file($folderName."/".$fileToDl);
          die "ERROR [$?]: Annocript can't download $fileToDl. Please try to download the file manually in the $folderName folder using:".
          "'wget  $folderName."/".$fileToDl' and restart Annocript using the same session. Sorry for this drawback. ".$res->status_line, "\n"; 
        }
      }
  }
}

=head2 compare_dimensions 

 Title   : compare_dimensions
 Usage   : compare_dimensions( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.
 
 Returns : 1 if succeed, -1 otherwise

=cut
sub compare_dimensions{
  my $originalDim = shift;
  my $filePath = shift;
  
  my $retVal = -1;
  
  if ( (-s $filePath) == ($originalDim/1024) ) {
      $retVal = 1;
  }
  
  return $retVal;
}


=head2 extract_file
 Title   : extract_file
 Usage   : extract_file( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen

 Returns : nothing
=cut
sub extract_file{
  my $input = shift;
  my $outDir = shift;
	
	my $ae = Archive::Extract->new( archive => $input);
  $Archive::Extract::PREFER_BIN = 1;
  
  my $resp = $ae->extract( to => $outDir );
  my $path = File::Spec->catfile( $outDir, $ae->files->[0] );
  
  #If there are some problems delete the file and die
  if ( !$resp ){
     delete_file($path);
     die "ERROR [$?]: ".$ae->error;
  } 
  
  #my $path = File::Spec->catfile( $outDir, $ae->files->[0] );
  printf "Decompression of ".$input." finished\n";
  
  return $path;
}

=head2 my_extract_file
 Title   : my_extract_file
 Usage   : my_extract_file( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen

 Returns : nothing
=cut
sub my_extract_file{
  my $input = shift;
  my $outDir = shift;
  
  my $command = '';
  my $outName = '';
  
  my $removeText = 'Please go in the folder and remove all the files coming out from this erroneous execution and try again.'.
  'If you get again an error please consider to unzip manually the files and leave them inside the folder. Then restart Annocript.';
  if ($input =~ /\.tar.gz$/){
      $command = "tar -zxf $input -C $outDir";
      try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText ?\n" ;
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText?\n" ;
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText?\n" ;
      }else{
        extract_file($input,$outDir);
      } 
      
  printf "Decompression of ".$input." finished\n";
}


=head2 extract_name

 Title   : extract_name
 Usage   : extract_name( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: extract the name from a complete path of a file. Even the file name only with the extension
              0: the complete name with the extension
              1: the name only
              2: the first two words joined by a dot
              gz: the name from the .gz
              targz: the name from the .tar.gz
              zip: the name from the .zip
              tar: the name from the .tar
              

 Returns : the name only

=cut
sub extract_name {
  my $filePath = shift;#Path to the file
  my $type = shift;#The type of file
  
  #Separate the path in pieces using slashes
  my @list = split("/",$filePath);
  my $complName = pop(@list);
  
  
  my @nameElements = split (/\./, $complName);
  my $name;
  switch ($type) {
    case ["0"]{ $name = $complName;}
    case ["1"]{ $name = $nameElements[0];}
		case ["2"]{ $name = $nameElements[0].'.'.$nameElements[1];}
    case ["gz"]{ $complName =~ /(\S+).gz/;
                 $name= $1;##Name for uncompressed file
                 }
		case ["targz"]{$complName =~ /(\S+).tar.gz/;
                 $name= $1;##Name for uncompressed file
                 }
    case ["zip"] {$complName =~ /(\S+).zip/;
                 $name= $1;##Name for uncompressed file
                 }
		case ["tar"] {$complName =~ /(\S+).tar/;
                 $name= $1;##Name for uncompressed file
                 }
    else { die	"ERROR [$?]: $type is not a valid input extracting a name: ?\n";}	
  }    
  return $name;
  
}



=head2 append_file_2_file

 Title   : append_file_2_file
 Usage   : append_file_2_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append the file on the path1 to the one on the path2
 
 Returns : nothing

=cut
sub append_file_2_file{
  my $path1 = shift;
  my $path2 = shift;
  
  sleep 1;
  open(DATA1, "<$path1") or annoDie( "Couldn't open file $path1");
  open(DATA2,">>$path2") or annoDie( "Couldn't open file $path2");
  
  while(<DATA1>){
   print DATA2 "$_";
  }
  close(DATA1);
  close(DATA2);
}

=head2 print_file

 Title   : print_file
 Usage   : print_file( - path1 -> the file to be print
                        )

 Function: will print the file on the screen
 
 Returns : nothing

=cut
sub print_file{
  my $path1 = shift;
  
  open(DATA1, "<$path1") or annoDie( "Couldn't open file $path1");
    
  while(<DATA1>){
   print "$_";
  }
  close(DATA1);
}


=head2 check_presence
 Title  : check_presence
 Usage  : check_presence(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : nothing

=cut
sub check_presence{
 my $fileToCheck = shift;
 
 my $retVal = 1;
 #print $fileToCheck." \n";
 if ( $fileToCheck ne ''){
   if(-z "$fileToCheck"){
    $retVal = 0;
    print "$fileToCheck is empty...!\n";
   }
   if(!-e "$fileToCheck"){
    $retVal = 0;
    print " $fileToCheck don't exists..\n";
   }
 }else{die ("\nERROR [$?]: file $fileToCheck has not been created: ?\n");}
 
 return $retVal;
}


=head2 check_presence
 Title  : check_presence
 Usage  : check_presence(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : nothing

=cut
sub check_presenceOld{
 my $fileToCheck = shift;
 
 my $retVal = 1;
 
 if ( $fileToCheck ne ''){
   if(-z "$fileToCheck"){
    die ("\nERROR [$?]: $fileToCheck is empty: ?\n");
   }
   if(!-e "$fileToCheck"){
    die ("\nERROR [$?]: $fileToCheck does not exist: ? \n");
   }
 }else{die ("\nERROR [$?]: file $fileToCheck has not been created: ?\n");}
 
 return $retVal;
}



=head2 join_files_with_cat

 Title   : join_files_with_cat
 Usage   : join_files_with_cat( -file1 = first file
                              - $file2 = second file
                              - outFile = otuput file
                               );

 Function:  Join two files in input and write them on a third file using the system call CAT

 Returns : nothing

=cut
sub join_files_with_cat{
 my $file1 = shift;
 my $file2 = shift;
 my $outFile = shift;
  
 #If the file is already present and its dimensions equals the sum of the two file it will not create again  
 unless ( (-e $outFile) and ( (-s $outFile) == (-s $file1) + (-s $file2) ) ){
    print "Concatenating ".$file1." and ".$file2." in ".$outFile."...";
   my $command="cat ".$file1." ".$file2." > ".$outFile;
    
   try_exec_command( qq($command) )== 1
     or die "ERROR [$?]: The concatenation of $file1 and $file2 failed: ?\n";
     
   print " ..completed. \n";
  }else{print "$outFile already present. Annocript will not create..\n";} 

 #close($outFile);
}


=head2 extract_special_chars

 Title  : extract_special_chars
 Usage  : extract_special_chars( - filePath => 'the path to the file');

 Function: this function searches special chars inside a file and returns a string with them
  
  Returns : a string

=cut		
sub extract_special_chars {
  my $filePath = shift;
  my $goodLetters = shift;
      
  my @collected = ();  
  open(FILE,"<$filePath");
        
  while(<FILE>){
    my $line = $_;
    $line =~ s/[$goodLetters]//;
    $line =~ s/\s//; 
  } 
  close(FILE);
} 
################################SPECIFIC FUNCTIONS FOR ANNOCRIPT MODULES


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

################################DATABASE FUNCTIONS


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

=head2 db_present

 Title   : db_present
 Usage   : db_present( -dbToSearch => 'name of the database);

 Function:  checks if a database exists.

 Returns : nothing

=cut
sub db_present {
  my $dbToSearch = shift;
  my $platform = shift;
	my $host = shift;
	my $port = shift;
	my $mySqlUser = shift;
	my $mySqlPass = shift;
  
  my $dsn = "dbi:$platform:information_schema:$host:$port";
              
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1} );
  my $databases = $dbh->selectcol_arrayref("show databases like '".$dbToSearch."'");
  
  my $numKeys = @$databases;
  
  return $numKeys;
}



####################################FUNCTIONS TO OPERATE ON FASTA FILES



=head2 detect_fasta_type

 Title   : detect_fasta_type
 Usage   : detect_fasta_type( - configFilePath = path of the config file
                               );

 Function:  this subroutine controls the first line of the first sequence of a fasta file.
            If it is composed of only ATGCN then they are nucleotidic sequences.
            
 Returns :  A string either 'dna','rna','protein'.

=cut
sub detect_fasta_type{  
  my $fastaFilePath = shift;

  my $seqio  = Bio::SeqIO->new( -file => $fastaFilePath);
   
  my $obj = $seqio->next_seq();
   
  return $obj->alphabet;

}


=head2 count_sequences
 Title  : count_sequences
 Usage  : count_sequences( - fasta => 'the path to fasta file of which we want to count seqs' );

 Function: 	counts the sequences in a fasta file. Using grep it searches for the '>' on each line.
            at the first occurence on the line it goes on the successive.

 Returns : nothing

=cut
sub count_sequences {
  my $fasta = shift;
  
  my $totSeqs = `grep -c '>' $fasta`;
   
  return $totSeqs; 
}



=head2 create_slices

 Title  : create_slices
 Usage  : create_slices( - fasta => 'name of the fasta files',
                         - sliceName => 'name of the slice (also to create a folder)',
                         - numSlices => 'number of slices to be used. This will be also the number of FASTA file',
                         - sessionFolder => 'the folder where to put these files'
                      );

 Function: counts the total number of bases in the fasta file and creates $numSlices FASTA files with approximately the same number of
            bases.
			
 Returns : nothing

=cut
sub create_slices{
  my $fasta = shift;
  my $sliceName = shift;
  my $numSlices = shift;
  my $sessionFolder = shift;
  
  #counts the total number of bases in the fasta file 
  my $total_bases = calculate_total_bases($fasta);

  #In each slide there will be the same number of bases
  my $bases_per_slice = int($total_bases/$numSlices)-1;


  print "\nTotal bases in $fasta: $total_bases\n";
  print "The average number of bases in the $numSlices slices will be: $bases_per_slice\n\n";

  #The original fasta file
  my $seqin = Bio::SeqIO->new(-file => $fasta,
                              -format => 'fasta');

  my $slicesDir = $sessionFolder.'/'.$sliceName;
  
  #print "I'll create a new folder: $slicesDir\n";
  #creates a folder with the fasta file name
  unless(-d $slicesDir){
    print $slicesDir." doesn't exists. Creating folder...\n";
    mkdir $slicesDir or die "ERROR [$?]: can't create folder $slicesDir: ?\n ";
  }
  
  die "Unable to change directory in $slicesDir.\n" unless chdir $slicesDir ;
  
  my $written = 0;#Number of bases written
  my $c = 1;#Number of files written

  my $name = "$sliceName\_$c.fasta";
  my $seqout;
    
  #Creates a SeqIO object for a new fasta file
  $seqout = Bio::SeqIO->new(-file => ">$name",
                            -format => 'fasta');

  #for each sequence in the original fasta
  while(my $seq = $seqin->next_seq) {
   # if the bases written are more then the decided threshold.. stop the object construction
   #and creates another one
    reset_seqout(\$c,\$written,\$seqout,$sliceName) if $written > $bases_per_slice;
    $seqout->write_seq($seq);
    $written += $seq->length;
  }
}


=head2 calculate_total_bases

 Title  : calculate_total_bases
 Usage  : calculate_total_bases( - fasta => 'name of the fasta file',
                      );

 Function: This subroutine counts the total number of bases in a fasta file by summing the lengths of each sequence.
            It is called

 Returns : nothing

=cut
sub calculate_total_bases {
  my $fasta = shift;
  
  my $total_bases = 0;
  my $seqin = Bio::SeqIO->new(-file => $fasta,
                              -format => 'fasta');
  
  while(my $seq = $seqin->next_seq) {
   $total_bases += $seq->length;
  }
  $seqin->close;
  
  return $total_bases;
}


=head2 reset_seqout

 Title  : reset_seqout
 Usage  : reset_seqout( - command => 'the command to execute',
                      );

 Function:  if the bases written are more then the decided threshold.. stop the object construction
              and creates another one
 Returns : nothing

=cut
sub reset_seqout {
  my ($c) = shift;
  my ($written) = shift;
  my ($seqout) = shift;
  my $sliceName = shift;
  
  $$c++;#increase counter for files
  $$written = 0;#initialize the number of bases written
  $$seqout->close;#Close the object

  #Starts a new file
  my $name = "$sliceName\_$$c.fasta";
  $$seqout = Bio::SeqIO->new(-file => ">$name",
                            -format => 'fasta');
}

########################################SYSTEM OPERATIONS


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

=head2 try_exec_command

 Title   : try_exec_command
 Usage   : try_exec_command( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.
 
 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_commandOld{
    my $command = shift;
    
    my $maxTimes = 5;
    my $success = -1;
    my $timesCount = 0;
    
    while ($success == -1 and $timesCount < $maxTimes){
        if ( (system $command) == 0) {
          $success = 1;
        }
        else{
         print "ERROR [$?]: an error occurred while executing $command. Trying again... $?\n"; 
         $timesCount++;
        }
    }
    return $success;
}


=head2 num_processors

 Title   : num_processors
 Usage   : num_processors(  - array = the array to print
                               );

 Function:  gives the number of processors of the system. It reads on the file cpu info that has written inside
            'processor' for each processor.
            
 Returns : nothing

=cut
sub num_processors {
  open FILE, "/proc/cpuinfo";
  my $num_cpus = 0;
    while(<FILE>) {
        if (/processor/) { $num_cpus++; }
    }
    close FILE;

  return $num_cpus;
}

=head2 ram_memory

 Title   : ram_memory
 Usage   : ram_memory(  - array = the array to print
                               );

 Function:  gives the ram memory of the system
 
 Returns : amount of ram in Bytes

=cut
sub ram_memory {
  open FILE, "/proc/meminfo";
  my $ram = "";
    while(my $row = <FILE>) {
      my $string = "MemTotal:"; 
      
        if ($row =~ /$string/) { 
          $row =~ s/$string//g;
          $row =~ s/kB//g;
          $row =~ s/\s//g; 
          $ram = $row;
          last;
        }
    }
    close FILE;
 #print $ram;#DEBUGCODE
 return $ram;
}


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

=head2 is_folder_empty

 Title  : is_folder_empty
 Usage  : is_folder_empty( - dirname => 'the folder to check',
                      );

 Function: Checks if a folder is empty

 Returns : nothing

=cut
sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

#################################CONTROLS ON STRINGS OF THE CONFIGURATION FILES

=head2 checkLink

 Title   : checkLink
 Usage   : checkLink( - configHash = the piece of the config hash for DB_CREATION
							 - local folder = the folder where actually is the DB_CREATION script
                               );

 Function:   check  if a given string in input is a correct link. Otherwise advice the user of the fact that it is not a correct link.
 Returns : nothing

=cut
sub checkLink {
    my $link = shift;
    my $retVal = 0;
    
    if ($link =~ /[http|shttp|ftp|sftp]:\/\/.+\/.*/){
        $retVal = 1;
    } 
    return $retVal;
 } 




=head2 correct_type
Title  : correct_type
 Usage  : correct_type( -number => 'a number to check',
                      -typeWanted => 'the expected type');

 Function: 	 Check if the number is of the type given in input
  
  Returns 1 if true

=cut
sub correct_type {
  my $number = shift;
	my $typeWanted = shift;
	
  my $type = "";
	my $ret = 0;

  if ($typeWanted eq "positiveint"){
    #contains onli digits
    if ( $number=~ m/^\d+$/) {
      $ret=1;
    }
  }elsif ($typeWanted eq "real"){
    #Contains digits starting with (+-), separated by '.' and can be added the 'E-x'
    if ( $number =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ){
        $ret=1;
    }
  }
  return $ret;
}

###########################OPERATIONS ON DATA STRUCTURES

=head2 deallocate_hash

 Title   : deallocate_hash
 Usage   : deallocate_hash(  - hash = the reference to the hash to deallocate (\@hash)
                               );

 Function: go through all the hash and deallocates the structure
 Returns : nothing

=cut
sub deallocate_hash{
  my ($hash) = shift;
  
  #Try to deallocate the huge hash
  foreach my $key (%$$hash){ 
    undef $key
  }
  undef $$hash;
}



=head2 print_array

 Title   : print_array
 Usage   : print_array(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints simply the values inside the array in input. Give it as \@array. USes the join function.
 Returns : nothing

=cut
sub print_array{
  my $array = shift;
  
  
  print join ' ', @$array, "\n";
  
  print "\n";
}


=head2 save_hash
 Title   : save_hash
 Usage   : save_hash(  - outHashPath -> the complete path of the file where the hash has to be written
                      - hashToWrite -> the hash to write
                      - hashName -> the hash name to use);

 Function:  this subroutine saves an hash given in input into a file
 
 Returns : nothing

=cut
sub save_hash {
  my $outHashPath = shift;
  my $hashToWrite = shift;
  my $hashName = shift;
          
  if (! open(BXHASH,">$outHashPath")){ return ("ERROR: Failure opening '$outHashPath' - $!");}
    #annoPrint (Dumper($hashToWrite); #DEBUGCODE
    #annoPrint (Data::Dumper->Dump([$hashToWrite], [$hashName]);#DEBUGCODE
    print BXHASH Data::Dumper->Dump([$hashToWrite], [$hashName]);
    close(BXHASH);
}
  
  
########################################INTERNET OPERATIONS

=head2 ftp_fetch_file

 Title   : ftp_fetch_file
 Usage   : ftp_fetch_file( -fileAdd => is the path of the file
                          -folderName => is the folder where we have to download
                          -user => username froruniprot
                            -pass => password for uniprot
                               );

 Function:  Given in input a path to a file, it download it in a specific folder. 

 Returns : nothing

=cut
sub ftp_fetch_file{
  my $fileAdd = shift;
  #print "file name: ".$fileAdd."\n";
  my $folderName = shift;
  my $user = shift;
  my $pass =shift;

  my $result = 1;
  
  my $maxTimes = 5;#MAx number of times to execute the command
  my $success = -1;#If it is success finish trying
  my $timesCount = 0;#Counter of faults
   
  #attempts maxTimes to access and download the file  
  while ( ($success == -1) and ($timesCount < $maxTimes) ){
   ## build a File::Fetch object ###
   my $ff = File::Fetch->new(uri => $fileAdd);
      # fetch the uri to folderName ###
      if ( $ff->fetch(to =>$folderName) ){
        $success = 1;
      }
      else{
        print "ERROR [$?]: an error occurred while fetching file $fileAdd. Trying again...\n"; 
        $timesCount++;
      }
  }
  #print $ff->error();
  return $success;
}

=head2 my_head

 Title   : my_head
 Usage   : my_head( - link => the link to check
                    );

 Function: It checks if a link to external resource is working. If it is it returns 1. If it is not then it will check
          what kind of error is coming out and alert the user with a message.

 Returns : 1 if link is good, 0 otherwise

=cut
sub my_head{
    my $link = shift;
    my $email = shift;
    
    my $retVal = 0;
    
    print "Checking $link\n";
    #If the link is uniprot then does a particular control
    #if ($link =~ /uniprot.org/){
      #$retVal = uniprot_access($link,$email);
    #}else{
      my $req = HTTP::Request->new(GET => $link);
      my $UA = new LWP::UserAgent;
      my $res = $UA->request($req);
       
      if($res->is_success) {
          #print $res->as_string();
          #$retVal = $res->header('Content-Length');
          $retVal = 1;
      } else {
          if ( (split(/ /,$res->status_line))[0] eq '404'){
            print "WARNING:(404 File doesn't exists). If you are sure of this link leave Annocript running.\n";
          }
          #print "Error: " . $res->status_line . "\n";
          $retVal = -1;
      }      
    #}
  return $retVal;
}


=head2 uniprot_access

 Title   : uniprot_access
 Usage   : uniprot_access( - link => the link to check
                    );

 Function: THIS FUNCTION I WAS USING WHEN IN THE INSTITUTE I WAS HAVING THIS STRANGE ERROR 404 ONLY FOR UNIPROT LINKS.
          BUT IF IT IS USELESS IT CAN BE ERASED
      It checks if a Uniprot link to external resource is working. If it is it returns 1. If it is not then it will check
          what kind of error is coming out and alert the user with a message.

 Returns : 1 if link is good, 0 otherwise

=cut
sub uniprot_access{
    my $link = shift;
    my $contact = shift;
    
    my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
    my $res = $agent->get($link);

    my $retVal = 0;
    
    if($res->is_success) {
        #print $res->as_string();
        $retVal = 1;
    } else {
        if ( (split(/ /,$res->status_line))[0] eq '404'){
          print "Uniprot.org gave strange response for: $link. Annocript will continue and try to download the file the same.\n";
        }
        #print "Error: " . $res->status_line . "\n";
        $retVal = -1;
    }
    return $retVal;
}
