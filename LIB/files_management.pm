#Copyright (C) <2018>  <Francesco Musacchia>

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
    
package LIB::files_management;
## files_management.pm
#Author: Francesco Musacchia  (2018)
#Permits the management of an files with many manipulation functions
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( read_and_print check_presence 
								compare_dimensions delete_directory 
								dl_and_extract download_file extract_archive
								extract_special_chars join_files_with_cat match_array_with_string
								my_extract_any_file my_extract_file_sys print_file check_presence
								extract_file_folder merge_columns extract_col_from_file 
								create_folder is_folder_empty file_list_to_array
								search_file_with_pattern file_not_present
								join_files_with_zcat save_hash load_hash extract_columns_from_file
								append_str_2_file file_num_rows invert_cols_position get_col_index
								delete_file insert_col_in_file_table delete_columns
								extract_colnum_from_file_linux file_name_wostrange_chars
								append_hash_to_file append_string_at_line_to_file 
								substitute_start_of_file check_URL delete_rows_containing
								shift_column compress_folder list_to_array insert_line_into_file
								separate_elements_on_col tar_folder_sys
								checkOutputsPresence count_lines_file extract_name);
}
use strict;
use warnings;
#FILE MANAGEMENT MODULES
use File::Path qw(make_path remove_tree ) ;
use File::Fetch; #TO use FTP download functions
use File::Copy; #Used to move files
use List::Util qw(first);#to get the index of a string in an array

use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;
use Archive::Extract;
use Archive::Tar;
#use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use IO::Compress::Gzip qw(gzip $GzipError);
use Data::Dumper;#To print the hashes
#use LIB::utilities qw(  extract_fasta_from_fasta );
use Storable;#To store data structures


#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(try_exec_command print_and_log log_and_exit);

#Using a library for standard utilities								
use LIB::std_lib qw(correct_type);


#################GENERIC OPERATIONS ON FILES


=head2 extract_name

 Title   : extract_name
 Usage   : extract_name( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: extract the name from a complete path of a file. Even the file name only with the extension
              0: the complete name with the extension
              1: the name only
              2: the first two words joined by a dot
              noext: just remove the last extension from the path
              no2ext: remove the last 2 extensions from the path
              gz: the name from the .gz
              targz: the name from the .tar.gz
              zip: the name from the .zip
              tar: the name from the .tar
              fqgz: the name before two extensions

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
  if ($type eq "0"){ $name = $complName;}
  elsif ($type eq "1"){ $name = $nameElements[0];}
  elsif ($type eq "2"){ $name = $nameElements[0].'.'.$nameElements[1];}
  elsif ($type eq 'noext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq 'no2ext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq "gz"){ $complName =~ /(\S+).gz/;
                 $name= $1;##Name for uncompressed file
                 }
  elsif ($type eq "targz"){$complName =~ /(\S+).tar.gz/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "zip"){$complName =~ /(\S+).zip/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "tar"){$complName =~ /(\S+).tar/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "fqgz"){
		$complName =~ /(\S+)\.\S+\.\S+/;
		
                 $name= $1;##Name for uncompressed file
                 }
   else { die	"ERROR [$?]: $type is not a valid input extracting a name: ?\n";}
  return $name;
  
}		
		
=head2 shift_column

 Title  : shift_column
 Usage  : shift_column ();
 
 Function: Shifts the column "col2shift" after the column "afterthis" into the file "intable".
					col2shift and afterthis must be the strings 
  
 Returns : the final output

=cut		
sub shift_column {
	my $intable = shift;
	my $col2shift = shift;
	my $afterthis = shift;
	my $log_file = shift;
	
	#Get the position of the ivf field to remove it later
	my $col2shift_ind = get_col_index($intable,$col2shift);	
	#Set the file names for the column which will be extracted from the file
	my $col_f = $intable.".temp_col2shift";
	
	#Extract the columns into the file
	print_and_log( "Extracting column $col2shift from $intable to $col_f\n",$log_file);#DEBUGCODE
	extract_col_from_file($intable,$col2shift,$col_f);
	
	#Set the array that will be the column to insert and insert it
	my @list = list_to_array($col_f,'NO_NEW_LINE');			
	#	push the header into the list
	unshift @list, $col2shift;
	
	#Choose the point where to insert and insert it

	my $afterthis_ind = get_col_index($intable,$afterthis);
	print_and_log( "Inserting column $col2shift from file $col_f into $intable after $afterthis [$afterthis_ind] \n",$log_file);
	print_and_log( "list array has ".scalar(@list)." elements and the first is ".$list[0]."\n",$log_file);#DEBUGCODE
	insert_col_in_file_table($intable,\@list,$afterthis,'i');
	#Delete the previous columns using the previously obtained index $col2shift_ind
	#If you put them before the previous position, now that position needs you sum 1
	my @cols_to_remove = ();
	if ( $afterthis_ind < $col2shift_ind){
		@cols_to_remove = (($col2shift_ind + 1));
	}else{
		@cols_to_remove = (($col2shift_ind));
	}
	print_and_log( "Deleting ".scalar(@cols_to_remove)." [".$cols_to_remove[0]."] columns from $intable\n",$log_file);#DEBUGCODE
	delete_columns($intable,\@cols_to_remove);				
		
}	



=head2 file_name_wostrange_chars

 Title   : file_name_wostrange_chars
 Usage   : file_name_wostrange_chars(  );

 Function:  this subroutine check if a filename with a specific extension
					does not have dots inside and if contains only chars, numbers or _ and -
            Needs that the file in input has an extension!
            
 Returns : 1 if check is positive, 0 otherwise

=cut
sub file_name_wostrange_chars{
	my $file_name = shift;
	my $extension = shift;
	
	my $retVal = 1;
	
 #The name has to be without dots inside. Only one dot can separate the extension from the name
  my $dots=0;
  my @line = split (//,$file_name);
  foreach my $char(@line){
      if ($char eq '.'){
          $dots++;
      }
  }
  if ($dots>1){
    die "Please change the file name: $file_name by removing all internal dots (not the one used for the extension)...\n";
		$retVal = 0;
  }
  
 
 #print "File name: ".$file_name."\n";
  my @name_parts = split(/\./,$file_name);

  #It cannot be longer than 50 chars
  if (@name_parts > 1){
		my $name = $name_parts[0];
    #Checks for permitted characters
    if ( length($name) > 50 or ($file_name !~ /^[A-za-z0-9\_\-]+\.($extension)$/i) ){
       die "$name is not a correct file name. Allowed characters [A-Z,a-z,0-9,_,-]. Allowed extensions [.$extension]. Max length: 50 chars\n";
			$retVal = 0;
    }
	}
	return $retVal;
}
 
 
 
=head2 extract_colnum_from_file_linux

 Title  : extract_colnum_from_file_linux
 Usage  : extract_colnum_from_file_linux( - filePath => 'the path to the file');

 Function: given a file and a name of a column (or the index), fetches all the values
					under the selected column
  
  Returns : a file with only one column

=cut		
sub extract_colnum_from_file_linux {
	my $file = shift;
	my $col_num = shift;
	my $retFile = shift;
	my $sep = shift;
	
	my $command = "cut -f".$col_num." -d '".$sep."' $file > $retFile";
	system($command) == 0 or die "Unable to execute command: $command\n";
}


=head2 read_and_print

 Title   : read_and_print
 Usage   : read_and_print( -file => 'a path to a file to print'
			       );

 Function: Reads a file using the more command

 Returns : nothing

=cut
sub read_and_print {
		my $file = shift;
		my $command = "more $file";
		try_exec_command($command) or die "Unable to execute command: $command\n";
}

=head2 file_list_to_array

 Title   : file_list_to_array
 Usage   : file_list_to_array(    );

 Function:  puts lines of a file inside an array and returns the array
 Returns : nothing

=cut
sub file_list_to_array{ 
	my $listFile = shift;
	
	my @list = ();
	
	open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
	while ( my $row = <LIST>){
			chomp($row);
			if (-e $row ){
				push(@list,$row);
			}else{
				die "ERROR in $listFile: $row does not exist.\n ";
			}
	}
	close(LIST);
	return @list;
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

=head2 search_file_with_pattern

 Title   : search_file_with_pattern
 Usage   : search_file_with_pattern( -pattern => 'a pattern to search'
													-folder => the folder where the file with pattern is to search
			       );

 Function:  Opens a directory and searches a filename containing the pattern
				given in input. Gives an error if more than one file is found. In that
				case the first found is given
 Returns : the file name identified or an empty string

=cut
sub search_file_with_pattern {
 my $pattern = shift;
 my $folder = shift;
 
 opendir DIR, $folder or die "ERROR [$!]: cannot open dir $folder\n";
 my @files= readdir DIR;
 closedir DIR;
 
 #Search in folder
 my $ret_file = "";
 my @results = grep {/$pattern/} @files;
  if ( scalar(@results) == 1 ){
		$ret_file = $results[0];
	}elsif ( scalar(@results) > 1 ){
			print "WARNING: too many files matching $pattern in $folder\n";
		}else{
			print "Can't find a file matching $pattern in $folder..\n";
			}
			
 return $ret_file;			
}


=head2 delete_file

 Title   : delete_file
 Usage   : delete_file( -filePath => 'a path to a file to delete'
			       );

 Function:  Delete a file in a given location. It can erase also more 
						than one file starting with the same name. This
            can be done by using the usual '*'. When it is found the 
            'glob' keyword will be used.
 Returns : Error code: 1 -> correctly deleted; -1 -> error with perl function; -2 file does not exist

=cut
sub delete_file{
  my $filePath = shift;
 
  my $retVal = 0;
  
	#With this first IF we can delete a set of file with the same name
	if ($filePath =~ /\*/){
		$retVal = 1;
		unlink glob $filePath or $retVal = -1;
	 }elsif ( -e $filePath ){
			if  ( unlink($filePath) == 1) {
				$retVal = 1;
				#deleted successfully
			}else{
			# not deleted for some problems with unlink subroutine or the file does not exist anymore
				$retVal = -1;
			}
	}else{
	# does not exist
		$retVal = -2;
	}
    
  return $retVal;
}

=head2 delete_directory

 Title   : delete_directory
 Usage   : delete_directory( -filePath => 'a path to a directory to delete'
			       );perl  VarGenius/get_sample_sheet.pl --help

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
  my $originalFileDim = shift;
  
  #Download the file 
  unless (-e $uncompressed){#Doesn't exists the  file, proceed
		unless (-e $compressed) {#Doesn't exists the compressed, download it
			 #Download
       #print "\n Downloading File: ".$link."\n";
        if (defined ($user) and defined ($pass)){
           download_file($link,$dbDataFolder,$user,$pass);
        }else{
           download_file($link,$dbDataFolder);
        }
        #Check if the dimensions are  equal to those taken during the link checking
        #If original is 0 it means that the check has been incorrect
        if ( defined $originalFileDim ){
					if ( $originalFileDim > 0){
						if (compare_dimensions($originalFileDim,$compressed) <= 0){
								die "Error: ".$compressed." has not been well downloaded. Delete it before to start again Annocript.\n";
								#delete_file($dbDataFolder."/".$compressed);
						}
					}
				}
		 }else{ print "File ".$compressed." already downloaded...\n"; }
		#Uncompression
    #print "\n Uncompress File: ".$dbDataFolder."/".$compressed."-->".$dbDataFolder."/".$uncompressed."\n";
    my_extract_any_file($compressed,$dbDataFolder);
 }else{print "File ".$uncompressed." already present...\n";}#Exists the file 
}



=head2 download_file

 Title   : download_file
 Usage   : download_file( -fileAdd => is the path of the file
                          -folderName => is the folder where we have to download
                          -user => username froruniprot
                            -pass => password for uniprot
                               );

 Function:  Given in input a path to a file, it download it in the data folder. It assigns as name of the file the last
#element of the array composed by the split of the string path.
The function tries more than one download types. The first one is with the GET of LWP. If it does not work, Annocript tries again with LWP Get
and if this again doesn't then it uses the Fetch::File class that calls various methods.
If this kind of download still does not work and we are trying to download Uniprot databases, then other domains are used because it can happen
that one of them is not working. The access is tried on Uniprot.org, EBI.org, Expasy.org

 Returns : it gives 1 if download succed, 0 if it is not, 2 if succeed with a change of link. In this case the new_link is valued

=cut
sub download_file{
  my $fileAdd = shift;
  #print "file name: ".$fileAdd."\n";
  my $folderName = shift;
  my $user = shift;
  my $pass =shift;
  my $originalFileDim = shift;
  
  my $new_link = '';
  my $retVal = 0;

  print "Downloading: $fileAdd...";
  
  my $ua = LWP::UserAgent->new;
  my $req = HTTP::Request->new (GET => $fileAdd);
  
  if (defined($user) and defined ($pass)){
    $req->authorization_basic($user, $pass);
  }
  my $fileToDl = extract_name($fileAdd,"0");
  my $res = $ua->request($req, $folderName."/".$fileToDl);

  if ($res->is_success) {
     print "..completed!\n";
     $retVal = 1; 
  }
  else {
	  #Here we add a specific control and re-download because sometimes for problems of networks
	  #the file is not downloadable but it becomes after some seconds.
     #print "Some problem occurred while downloading $fileAdd. Trying again...\n";
	   my $res2 = $ua->request($req, $folderName."/".$fileToDl);
	   if ($res2->is_success) {
			print "..completed!\n";
      $retVal = 1;
      }
      #Commented because I do not need to download at this time. Uncomment if needed
      #else{
     
        ##This command will force the download with different methods
        #if (ftp_fetch_file($fileAdd,$folderName,$user,$pass) == 1){
            #print "..completed!\n";
            #$retVal = 1;
        ##If it does not work and we are searching Uniprot, then access different servers
        #}#elsif (match_array_with_string($fileAdd, \@uni_domains) == 1){
          ##$new_link  = download_from_diff_sources($fileAdd,$folderName,$user,$pass);
          ##$retVal = 2;
        ##}
      #}
  }
  if ($retVal == 0){
    die "Unable to download file $fileAdd. Restart again later. Annocript will close...\n";
  }
  return $retVal, $new_link;
}


=head2 match_array_with_string 

 Title   : match_array_with_string
 Usage   : match_array_with_string( -string => the string to search in
                                    - @array => the array with elements to search
                               );

 Function: Searches if at least one element of an array is present in a string
 
 Returns : 1 if succeed, 0 otherwise

=cut
sub match_array_with_string {
  my $string = shift;
  my $array = shift;

  my $found = 0;
  
  foreach my $el (@$array){
    if ( $string =~ /\Q$el\E/){
      $found = 1;
    }
  }
  return $found;
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

=head2 compress_folder
 Title   : compress_folder
 Usage   : compress_folder( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of a file and that of the final compressed
            and compresses it as a tar

 Returns : nothing
=cut
sub compress_folder{
	my $dir = shift;
	my $infile = shift;
	my $outfile = shift;
	
	# Create a new tar object:
	my $tar = Archive::Tar->new();
	
	#Go to the input dir path
	chdir $dir;
	
	#Compress
	$tar->create_archive( $outfile,COMPRESS_GZIP,glob $infile."/*");
}

#=head2 gzip_folder
 #Title   : gzip_folder
 #Usage   : gzip_folder( -input = the file compressed
                            #- output = the file uncompressed
                               #);

 #Function:  takes in input the path of a file and that of the final compressed
            #and compresses it as a tar

 #Returns : nothing
#=cut
#sub gzip_folder{
	#my $dir = shift;
	#my $infile = shift;
	#my $outfile = shift;
	
	#my $in = $dir."/".$infile;
	#my $out = $dir."/".$outfile;
	
	#my $zip = Archive::Zip->new();

  ## Add a directory
  #my $dir_member = $zip->addDirectory( $in);
 
  ## Save the Zip file
  #unless ( $zip->writeToFileNamed($out) == AZ_OK ) {
    #die 'write error';
  #} 
#}

=head2 tar_folder_sys
 Title   : tar_folder_sys
 Usage   : tar_folder_sys( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of a file and that of the final compressed
            and compresses it as a tar with the linux command

 Returns : nothing
=cut
sub tar_folder_sys{
	my $dir = shift;
	my $infold = shift;
	my $outfold = shift;
	
	
	#Use the system command TAR to tar the archive
	my $command = "tar -czvf $dir/$outfold $dir/$infold";
	#print_and_log( "Executing: $command..\n",$log_file);
	print "Executing: $command..\n";
	try_exec_command($command) or die "Unable to execute command: $command\n";	
}


=head2 extract_archive
 Title   : extract_archive
 Usage   : extract_archive( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen

 Returns : nothing
=cut
sub extract_archive{
  my $input = shift;
  my $outDir = shift;
	
  my $path = '';
	my $ae = Archive::Extract->new( archive => $input);
  $Archive::Extract::PREFER_BIN = 1;
  
  my $resp = $ae->extract( to => $outDir );
  $path = File::Spec->catfile( $outDir, $ae->files->[0] );
  
  #If there are some problems delete the file and die
  if ( !$resp ){
     delete_file($path);
     die "ERROR [$?]: ".$ae->error;
  } 
  
  #my $path = File::Spec->catfile( $outDir, $ae->files->[0] );
  #printf "Decompression of ".$input." finished\n";
  
  return $path;
}

=head2 my_extract_any_file
 Title   : my_extract_any_file
 Usage   : my_extract_any_file( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen.
            It can uncompress .tar.gz calling the extract_archive subroutint
            .gz, .bz2 and .zip calling anyuncompress subroutine.
            Thus, it uncompress always using PERL modules.

 Returns : nothing
=cut
sub my_extract_any_file{
  my $input = shift;
  my $outDir = shift;
  
  print "\n Uncompressing $input...\n";
  my $command = '';
  my $outName = '';
  if ($input =~ /\.tar.gz$/){
    extract_archive($input,$outDir);
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        #$command = "gunzip -c $input > $outName";
        anyuncompress $input => $outName
          or die " gz uncompress failed: $AnyUncompressError\n";
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        anyuncompress $input => $outName
          or die "zip uncompress failed: $AnyUncompressError\n";
      }elsif ($input =~ /\.bz2$/){
        $outName = extract_name($input,"bz2");
        $outName = $outDir."/".$outName;
        anyuncompress $input => $outName
          or die "bz2 uncompress failed: $AnyUncompressError\n";
        }  
  #print "Uncompression of ".$input." finished\n";
  print "...completed!\n";
}

=head2 my_extract_file_sys
 Title   : my_extract_file_sys
 Usage   : my_extract_file_sys( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen
                It extracts by using system calls.

 Returns : nothing
=cut
sub my_extract_file_sys{
  my $input = shift;
  my $outDir = shift;
  
  my $command = '';
  my $outName = '';
  
  my $removeText = 'Please go in the folder and remove all the files coming out from this erroneous execution and try again.'.
  'If you get again an error please consider to unzip manually the files and leave them inside the folder. Then restart Annocript.';
  if ($input =~ /\.tar.gz$/){
      $command = "tar -zxf $input -C $outDir";
      try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText \n" ;
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }else{
        extract_archive($input,$outDir);
      } 
      
  printf "Decompression of ".$input." finished\n";
}



=head2 extract_file_folder

 Title   : extract_file_folder
 Usage   : extract_file_folder( -filePath => 'complete path of the file',
			       );

 Function: extracts the folder of a file from its complete path
              

 Returns : the name only

=cut
sub extract_file_folder {
  my $filePath = shift;#Path to the file
  
  #Separate the path in pieces using slashes
  my @list = split("/",$filePath);
  pop(@list);
  my $folder = join("/",@list);
  return $folder;  
}



=head2 file_num_rows

 Title   : file_num_rows
 Usage   : file_num_rows( - path1 -> the file to be used)

 Function: returns number of rows of a file
 
 Returns : nothing

=cut
sub file_num_rows{
  my $path1 = shift;

 my $row_num = 0;
  open(DATA1, "<$path1") or die "Couldn't open file $path1";
    
  while(<DATA1>){
   $row_num++;
  }
  close(DATA1);
return $row_num;
}


=head2 invert_cols_position

 Title   : invert_cols_position
 Usage   : invert_cols_position( - file -> the file to write,
                                - colA => the column index A
                                - colB => the column index B 
                                )

 Function: Inverts the position of the column A and column B.
					If colA=10 and colB=11, The column 11 will be the column 10 and
					viceversa.
						 
 Returns : the same file is modified

=cut
sub invert_cols_position{
  my $file = shift;
  my $colA = shift;
  my $colB = shift;
  
  my $tempFile = $file."_TEMP";
  my $command = 'awk -F"\t" '."'".'BEGIN{OFS=FS;}{t=$i;$i=$j;$j=t;}1'."'"." i=$colA j=$colB $file > $tempFile";
	print "Executing command: $command\n";
	system ($command) == 0 or die "Unable to execute command: $command\n";
					
	move($tempFile,$file) or print "ERROR: unable to move $tempFile in $file\n";
	delete_file($tempFile);
}

  
=head2 append_str_2_file

 Title   : append_str_2_file
 Usage   : append_str_2_file( - path -> the file where to write,
                                - string => the string to be written)

 Function: 
						will append the string at the file path and put a newline.
					
 
 Returns : nothing

=cut
sub append_str_2_file{
  my $path = shift;
  my $string = shift;
  
  print "String: $string will be appended in $path..\n";
  
  #Check if the field exists
	open(DATA,">>$path") or die "Couldn't open file $path";
	print DATA $string."\n";
	close(DATA);

}

=head2 count_lines_file

 Title   : count_lines_file
 Usage   : count_lines_file( -filePath => 'a path to a file of which to count lines'
			       );

 Function:  Count lines in a file
 Returns : returns number of lines

=cut
sub count_lines_file {
	my $filePath =  shift;
	open(IN, "<$filePath");
	my @str = <IN>;
	close(IN);
	return scalar(@str);	
}

=head2 append_hash_to_file

 Title   : append_hash_to_file
 Usage   : append_hash_to_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append the hash in input on the file
 
 Returns : nothing

=cut
sub append_hash_to_file{
  my $file = shift;
  my $hash = shift;
  
  open(FILE, ">>$file") or die "Couldn't open file $file";
  
  foreach my $key (keys %$hash){
		print FILE $key." = ".$hash->{$key}."\n";
	}
  close(FILE);
}

=head2 append_string_at_line_to_file

 Title   : append_string_at_line_to_file
 Usage   : append_string_at_line_to_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: will append a string after a specific line (detected with a string)
					into a file. 
 
 Returns : nothing

=cut
sub append_string_at_line_to_file{
  my $file = shift;
  my $string_to_find = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	
	# print the lines before the change
	while( my $line = <FILE> ) {
    chop($line);
    print TEMP $line."\n";
    last if $line eq $string_to_find; # chosen string to find found
   }
   #Print the input line
	print TEMP $string_to_ins."\n";
	# print the rest of the lines
	while( my $line = <FILE> )  {
     print TEMP $line."\n";
   }
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
}

=head2 insert_line_into_file

 Title   : insert_line_into_file
 Usage   : insert_line_into_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: inserts an input line into a file given an index.
						The line is inserted at the index
					
 Returns : nothing

=cut
sub insert_line_into_file{
  my $file = shift;
  my $string_to_ins = shift;
  my $index = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	

	my $num_row = 0;
	# print the lines before the change
	while( my $line = <FILE> ) {
		if ( $num_row == $index){
			#Print the input string in the new file
			print TEMP $string_to_ins;		
		}
		#The current line is print after the new one
    print TEMP $line;
    $num_row++;
   }
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
	close(TEMP);
	close(FILE);
}

###SEPARATE the rows based on one column with multiple elemnts
##Example:
## chr1 9292 299923 GeneA,GeneB,GeneC
#becomes:
## chr1 9292 299923 GeneA
## chr1 9292 299923 GeneB
## chr1 9292 299923 GeneC
sub separate_elements_on_col{
	
	my $bed = shift;
	my $outname = shift;
	my $col2Sep_ind = shift;
	
	####################SEPARATE the rows based on one column with multiple elemnts
	open (FILE,"<$bed") or die "Cannot open $bed\n";
	open(NEWFILE,">$outname");
	my $numline = 0;
	#print "Executing gene separation on file $bed to fie $outname..\n";#DEBUGCODE
	foreach my $line (<FILE>){
		chop($line);
		my @fields = split("\t",$line);
		#print $line;#DEBUGCODE
		#Get the column to sep
		my $col2Sep = $fields[$col2Sep_ind - 1];
		my @elems = split(",",$col2Sep);	
		my $newline = "";
		foreach my $elem (@elems){
			#if (scalar(@elems) > 1){print "$elem \t";}
			for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
				$newline .= $fields[$col]."\t";
			}
			$newline .= $elem."\t";
			for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
				$newline .= $fields[$col]."\t";
			}			
			chop($newline);
			$newline .= "\n";
		}
		print NEWFILE $newline;
		#if (scalar(@elems) > 1){print "\n";}
		$numline++;
	}
	close(FILE);
	close(NEWFILE);	
}



#=head2 separate_elements_on_col

 #Title   : separate_elements_on_col
 #Usage   : separate_elements_on_col( - path1 -> the file to be copied,
                                #- path2 => the file to be written)

 #Function: given a file and an index of a column, separates the elements
						#separated with comma and writes a line for each element
					
 #Returns : nothing

#=cut
#sub separate_elements_on_col{
	#my $file = shift;
	#my $outname = shift;
	#my $col2Sep_ind = shift;
	
	#print  "Separating $file into $outname\n";
	#####################SEPARATE elements separated with comma
	#open (FILE,"<".$file) or die "Cannot open $file\n";
	#open(NEWFILE,">".$outname) or die "Cannot open $outname\n";
	#print  "Column $col2Sep_ind\n";
	
	#while (my $line = <FILE>){
		#print  "Column 2 $col2Sep_ind: $line\n";
		#chop($line);
		#my @fields = split("\t",$line);
		
		##If the field contains elements separated by comma
		#if ($fields[$col2Sep_ind] =~ /,/){
			
			##Get index of the column wanted
			#my @elements = split(",",$fields[$col2Sep_ind]);	
							
			#foreach my $element  (@elements){
				#my $newline = "";
				#for (my $col = 0; $col < ($col2Sep_ind - 1); $col++){
					#$newline .= $fields[$col]."\t";
				#}
				#$newline .= $element."\t";
				#for (my $col= $col2Sep_ind; $col < scalar(@fields); $col++){
					#$newline .= $fields[$col]."\t";
				#}			
				#chop($newline);
				#$newline .= "\n";
				#print NEWFILE $newline;
			#}			
		#}else{
				#print NEWFILE $line."\n";
		#}
	#}
	#print "Exiting\n";
	#close(FILE);
	#close(NEWFILE);	
#}


=head2 substitute_start_of_file

 Title   : substitute_start_of_file
 Usage   : substitute_start_of_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: substitutes the start of the file given in input.
					It searches for the string_to_find, and until that string
					writes only the string_to_ins in the new file. Once that the string
					is found it copies all the rest of the input file into the output
					
 Returns : nothing

=cut
sub substitute_start_of_file{
  my $file = shift;
  my $string_to_find = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	
	#Print the input string in the new file
	print TEMP $string_to_ins."\n";
	
	# print the lines before the change
	while( my $line = <FILE> ) {
    chop($line);
    last if $line eq $string_to_find; # chosen string to find found
   }

	print TEMP $string_to_find."\n";
	# print the rest of the lines
	while( my $line = <FILE> )  {
     print TEMP $line;
   }
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
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
  
  open(DATA1, "<$path1") or die "Couldn't open file $path1";
    
  while(<DATA1>){
   print "$_";
  }
  close(DATA1);
}


=head2 check_presence
 Title  : check_presence
 Usage  : check_presence(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : 1 if the file is ok, otherwise 0

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
    print " $fileToCheck does not exists..\n";
   }
 }else{
	 print ("The file to be checked is an empty string. Please check the code...\n");
	 $retVal = 0;
	}
 
 return $retVal;
}

=head2 file_not_present
 Title  : file_not_present
 Usage  : file_not_present(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : 1 if the filename is an empty string
					 2 if the file is not present
					 3 if the file is present but the size is 0

=cut
sub file_not_present{
 my $fileToCheck = shift;
 
 my $retVal = 0;
 #print $fileToCheck." \n";
 if ( $fileToCheck ne ''){
   if(-z $fileToCheck){
    $retVal = 1;
    print "$fileToCheck is present but empty...!\n";
   }
   #-e could return false if the directory is not readable by the program giving a misleading behaviour
	
   if(! (-e $fileToCheck) ){
    $retVal = 2;
		#print "$fileToCheck does not exists..\n";
   }
 }else{
	 print "The file to be checked is an empty string. Please check the code...\n";
	 $retVal = 3;
	}
 
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
    
   try_exec_command($command) > 0 
     or die "ERROR [$?]: The concatenation of $file1 and $file2 failed: ?\n";
     
   print " ..completed. \n";
  }else{print "$outFile already present. I will not create..\n";} 

 #close($outFile);
}



=head2 join_files_with_zcat

 Title   : join_files_with_zcat
 Usage   : join_files_with_zcat( -folder = folder with the files
                              - $out_zip = otuput file
                              - compr_f_list = list of compressed files
                               );

 Function:  Joins a list of compressed files in input in a single compressed
			file

 Returns : nothing

=cut
sub join_files_with_zcat{
	my $folder = shift;
	my $out_zip = shift;
	my $compr_f_list = shift;
	
	#$out_zip = $folder."/".$out_zip;
	if ( ! -e $out_zip ){
	my $temp_folder = extract_file_folder($out_zip);
		my $temp_out = $temp_folder."/temp.out";
		foreach my $comp_file (@$compr_f_list){
			#chomp($data_file);
			
			my $command =  "zcat ".$folder."/".$comp_file." >> $temp_out";
			print "Executing: $command\n";
			if (try_exec_command($command) < 1){
				die "Cannot execute $command\n";
			}
		}
		print "Now compressing with gzip to $out_zip \n";
		#Compress with gzip
		gzip $temp_out => $out_zip;
		#Remove the temp file
		delete_file($temp_out);
		print $temp_out." removed..\n";
	}else{
			die "ERROR: file $out_zip already exists. Please delete that before to start...\n";
	}
}

=head2 check_joined_with_zcat

 Title   : check_joined_with_zcat
 Usage   : check_joined_with_zcat( -folder = folder with the files
                              - $out_zip = otuput file
                              - compr_f_list = list of compressed files
                              - fusion_file = zipped file with the merge
                               );

 Function:  Joins a list of compressed files in input in a single compressed
			file

 Returns : nothing

=cut
sub check_joined_with_zcat{
	my $folder = shift;
	my $out_zip = shift;
	my $compr_f_list = shift;
	my $fusion_file = shift;
	
	
	$out_zip = $folder."/".$out_zip;
	if ( ! -e $out_zip ){
	
		my $temp_out = $folder."/temp.out";
		foreach my $comp_file (@$compr_f_list){
			#chomp($data_file);
			my $command =  "zcat ".$folder."/".$comp_file." >> $temp_out";
			print "Executing: $command\n";
			if (try_exec_command($command) < 1){
				die "Cannot execute $command\n";
			}
		}
		print "Now compressing with gzip to $out_zip \n";
		#Compress with gzip
		gzip $temp_out => $out_zip;
		#Remove the temp file
		delete_file($temp_out);
		print $temp_out." removed..\n";
	}else{
			die "ERROR: file $out_zip already exists. Please delete that before to start...\n";
	}
}






=head2 merge_columns

 Title  : merge_columns
 Usage  : merge_columns( - filePath => 'the path to the file',
												- outFile => 'output file',
												- separator => 'to indicate which char separates columns'
												);
 Function: all the columns are separated and merged row by row
  
  Returns : a string

=cut		
sub merge_columns {
	my $file = shift;
	my $outFile = shift;
	my $in_sep = shift;
	
	my $separator;
	
	print "Merging columns from $file...\n";
	if (!(defined $in_sep)){
		$separator = "\t";
	}else{
		$separator = $in_sep;
	}
	open (FILE,"<$file") or die "Cannot open $file\n";
	open (OUT,">$outFile") or die "Cannot open $outFile\n";
	
	while (my $row = <FILE>){
			my @pieces = split($separator,$row);

			print OUT join("\n",@pieces);
	}
	
	close(OUT);
	close(FILE);
}

=head2 extract_col_from_file

 Title  : extract_col_from_file
 Usage  : extract_col_from_file( - filePath => 'the path to the file');

 Function: given a file and a name of a column (or the index), fetches all the values
					under the selected column
  
  Returns : a file with only one column

=cut		
sub extract_col_from_file {
	my $file = shift;
	my $colName = shift;
	my $retFile = shift;
	my $sep_c = shift;
	
	open (COL_FILE, ">$retFile") or die ("Cannot open $retFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");
	
	my $header = <DATA>;
	
	my $sep = "\t"; 
	if (defined $sep_c){
		$sep = $sep_c;
	}	
	my @header_fields = split (/$sep/, $header);

	
	my $num_col = 0;
	
	#If the colName is the index of the column takes that
	if ( $colName =~ m/^\d+$/ ){
		$num_col = $colName;
	}
	#Otherwise it is the column name
	else{
		for my $field (@header_fields) {
			last if $field =~ m/$colName/;
			$num_col++;
		}		
	}
	
	while ( my $row = <DATA> ) {
		last unless $row =~ /\S/;
		chomp($row);
		my @fields = split /$sep/, $row;
		print COL_FILE $fields[$num_col]."\n";
	}

	close(DATA);
	close(COL_FILE);	
}




=head2 get_col_index

 Title  : get_col_index
 Usage  : get_col_index( - filePath => 'the path to the file');

 Function: given a file and a name of a column, searches the index
						of the column in input
  
  Returns : the index of the column as integer (starts from zero)
  
=cut		
sub get_col_index {
	my $file = shift;
	my $colName = shift;
	
	open (FILE, "<$file") or die ("Cannot open $file\n");
	
	my $header = <FILE>;
	
	my @header_fields = split (/\t/, $header);
	my $num_col = -1;
	
	for my $field (@header_fields) {
		$num_col++;
	  last if $field =~ m/$colName/;
	}
	
	close(FILE);	
	return $num_col;
}


=head2 extract_columns_from_file

 Title  : extract_columns_from_file
 Usage  : extract_columns_from_file( - filePath => 'the path to the file');
																			columns => a string with numbers separated by comma
 Function: given a file and a list of a columns indexes, returns only the 
					selected column. This script allows also to resort a table with the
					given order
					N.B. Starts from 0 as they are indexes of a PERL array
					
  Returns : a file with only the needed columns in the order as given in input

=cut		
sub extract_columns_from_file {
	my $file = shift;
	my $columns = shift;#a string with numbers separated by comma
	my $outFile = shift;
	
	open (COL_FILE, ">$outFile") or die ("Cannot open $outFile\n");
	open (DATA, "<$file") or die ("Cannot open $file\n");
	#Get the columns in an array
	my @cols = split(",",$columns );
	
	#Write a new file with the given columns only
	while ( my $row = <DATA> ) {
		#last unless $row =~ /\S/;
		chomp($row);
		my @fields = split /\t/, $row;
		if (scalar(@cols) > scalar(@fields) ){ die "ERROR in extract_columns_from_file: More columns ids than columns\n";}
		foreach my $col (@cols){
			print COL_FILE $fields[$col]."\t";
		}
		#Remove the last \t
		chomp($row);
		#print newline
		print COL_FILE "\n";
	}

	close(DATA);
	close(COL_FILE);	
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

=head2 create_folder

 Title  : create_folder
 Usage  : create_folder( - folderPath => 'the path to the folder');

 Function: creates the folder if it does not exist
  
  Returns : nothing
  
=cut		
sub create_folder {
	my $folderPath = shift;
	
	#Creation
	unless(-d $folderPath){
		print "$folderPath doesn't exists. Creating folder...\n";
		mkdir $folderPath or die "ERROR: can't create folder $folderPath\n ";
	}
}


=head2 insert_col_in_file_table
 Title   : insert_col_in_file_table
 Usage   : insert_col_in_file_table(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - an array of values in input
						# - a column name after which you want to add the array as the new column
						# opens the file, reads the header and finds the index in which should be 
						# inserted the column. Then it goes through all the rest of the file copying
						# the content in a new file and at that index it always insert the corresponding
						# value from the array set.
 
 Returns : nothing

=cut
sub insert_col_in_file_table {
	my $fileTable = shift;
	my ($newCol) = shift;#the vector to insert
	my $colNum = shift;#string representing the column
	my $job = shift;#Kind of job to exectue: substitution or insertion
	my $sep = shift;#separator to use, default: \t
	
	#If a separator is defined, use that, else it is tab
	unless (defined $sep) {$sep = "\t";}
 my $newFile = "newfile.temp";
 
 #Open input and output file
 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
 open (OUT,">$newFile") or die "Cannot open $newFile\n";
 
 my $row = <IN>;
 chomp($row);
 my @header = split($sep,$row );
 
 #This function will take the index of the first element in the header
 #which matches the string col
 my $index = first {$header[$_] eq $colNum} 0..$#header;
 
 print "$colNum is at position $index\n";
 my $totCols = scalar(@header);	
 my $row_ind = 0;
 
 #Here print the header
 if ($job eq 's'){
	print OUT join($sep,@header);
	print OUT "\n";
	
 }elsif ($job eq 'i'){
	 close (IN);
	 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
	 $index++;
	 }
 
 #For each line insert the given column
 while ( $row = <IN>){
		chomp($row);
		my @fields = split($sep,$row );
		#Take all the columns until the selected index
		my @toPrint = @fields[0 .. ($index-1)];
		#Print them separated with tab
		print OUT join($sep,@toPrint);
		#Then prints the column wanted
		print OUT $sep.$$newCol[$row_ind];#$newCol->[$row_ind]
		
		my $next_ind;
		if ($job eq 's'){
			$next_ind = $index+1;
		}elsif ($job eq 'i'){
			$next_ind = $index;
		}
		
		#Then finish to print the rest of the line
		if ($next_ind+1 <= $totCols) {
			print OUT $sep;
			@toPrint = @fields[$next_ind .. $#fields];
			print OUT join($sep,@toPrint);
		}
		print OUT "\n";
		
		$row_ind++;
	}
	
 close (IN);
 close (OUT);
 
	move($newFile,$fileTable) or print "ERROR: unable to move $newFile in $fileTable\n";
	delete_file($newFile);
}


=head2 delete_columns
 Title   : delete_columns
 Usage   : delete_columns(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - a set of column index that you want to remove (starting from 0)
						# opens the file, goes through all the file copying
						# the content in a new file and it will jump those indexes
 
 Returns : nothing

=cut
sub delete_columns {
	my $fileTable = shift;
	my $cols_to_rem = shift;#the vector to insert
	my $sep = shift;

	#If a separator is defined, use that, else it is tab
	unless (defined $sep) {$sep = "\t";}
 my $newFile = "newfile.temp";
 
 #Open input and output file
 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
 open (OUT,">$newFile") or die "Cannot open $newFile\n";
 

 #For each line insert the given column
 while ( my $row = <IN>){
		chomp($row);
		my @fields = split($sep,$row );
		my $col_num = 0;
		my $new_row = "";
		foreach my $field (@fields){
			if ( ! grep {/\b$col_num\b/} @$cols_to_rem){
				#Then prints the column wanted
				$new_row = $new_row.$field.$sep;
			}
			$col_num++;
		}
		chop($new_row);
		print OUT $new_row."\n";
	}
	
 close (IN);
 close (OUT);
 
	move($newFile,$fileTable) or print "ERROR: unable to move $newFile in $fileTable\n";
	delete_file($newFile);
}


=head2 delete_rows_containing
 Title   : delete_rows_containing
 Usage   : delete_rows_containing(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  #Given in input:
						#	- a file made as a table
						# - a symbol
						# opens the file, goes through all the file removing
						#all lines containing a given symbol in input
 
 Returns : nothing

=cut
sub delete_rows_containing {
	my $fileTable = shift;
	my $symbol = shift;#something to grep
	my $newFile = shift;
		
 #Open input and output file
 open (IN,"<$fileTable") or die "Cannot open $fileTable\n";
 open (OUT,">$newFile") or die "Cannot open $newFile\n";
 
 #For each line insert the given column
 while ( my $row = <IN>){
	 if ( $row !~ /$symbol/){
		print OUT $row."\n";		 
	 }
	}
 close (IN);
 close (OUT);

}


=head2 save_hash
 Title   : save_hash
 Usage   : save_hash(  - hash_file -> the complete path of the file where the hash has to be written
                      - hash -> the hash to write
                      );

 Function:  this subroutine saves an hash given in input into a file
					using the package Storable
 
 Returns : nothing

=cut
sub save_hash {
	my ($hash) = shift;
	my $hash_file = shift;
	#print "Storing hash..";
	#print $$hash->{'B'}."\n";#DEBUGCODE
	use Storable;
	store $$hash, $hash_file;
}

=head2 load_hash
 Title   : load_hash
 Usage   : load_hash(  - hash_file -> the complete path of the file where the hash is saved
                      );

 Function:  this subroutine loads an hash from a file
 
 Returns : the hash loaded

=cut
sub load_hash{
	my $hash_file = shift;
	return retrieve($hash_file);
}

=head2 check_URL

 Title   : check_URL
 Usage   : check_URL( - link => the link to check
                    );

 Function: It checks if an URL to external resource is working. If it is, it returns 1. If it is not then it will check
          what kind of error is coming out and alert the user with a message.

 Returns : 1 if link is good, -1 otherwise

    
=cut
sub check_URL{
    my $link = shift;
    
    my $agent = LWP::UserAgent->new();
    my $res = $agent->get($link);

    my $retVal = -1;
    
    if($res->is_success) {
        #print $res->as_string();
        $retVal = 1;
    } else {
        #if ( (split(/ /,$res->status_line))[0] eq '404'){
          #print "ERROR: ".$res->status_line.".\n";
        #}
        print "ERROR: ".$res->status_line."\n";
        $retVal = -1;
    }
    return $retVal;
}

=head2 list_to_array

 Title   : list_to_array
 Usage   : list_to_array( listFile = file path
													new_line = take or not the new line
										);

 Function:  puts lines of a file inside an array and returns the array.
					You can use it taking or not the new line character '\n' at the end of eac line
					by using the parameter new_line (NO_NEW_LINE, to take or nothing or NO to
					not take.
					DEFAULT: will take
 Returns : an array with the lines

=cut
sub list_to_array{ 
	my $listFile = shift;
	my $new_line = shift;
	
	my @list = ();
	
	#you can do it using the newline or not
	if ($new_line eq 'NO_NEW_LINE'){
		open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
		while ( my $row = <LIST>){
				chomp($row);
				push(@list,$row);
		}
		close(LIST);
	}else{
		#Reading all in one with <FILE> command
		open (FILE,"<$listFile") or die "Cannot open $listFile\n";
		@list = <FILE>;
		close(FILE);
	}
	return @list;
}

=head2 checkOutputsPresence

 Title : checkOutputsPresence
 Usage : checkOutputsPresence( -configHash => 'configHash',
					);

 Function: This function memorize if the output files are present and store the result in some variables of the config hash to keep it
          always known.

 Returns : YES if it is present, NO otherwise

=cut	 
sub checkOutputsPresence{
	my $outputName = shift;
	my $outFolder = shift;
	
	my $retVal = '';
	#Check the presence of the blastx output
	if (-e $outFolder.'/'.$outputName){
		 $retVal = 'YES';
	}else{ $retVal = 'NO';}
	return  $retVal;
}


1;
