
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
    
package LIB::html_management;
## html_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of HTML files for the pipeline
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( file_2_html html_page paths_hash_to_html var_stats_2_html
								table_2_html table_list_2_html testFTP
								ftp_fetch_file check_FTP_diff_sources
								checkLink dl_and_extract download_file
								check_url);
}
use strict;
use warnings;
use CGI; #load cgi module to produce HTML output
use File::Fetch; #TO use FTP download functions
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;
use Archive::Extract;

#HTTP FTP CHECKS
use LWP::UserAgent;
use LWP::Simple;
use HTTP::Request;
use HTTP::Response;
use HTTP::Date;
use Net::FTP;

#Using a library to manage files
use LIB::files_management qw(delete_file extract_name list_to_array my_extract_any_file
														match_array_with_string);


#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(print_and_log log_and_exit );

#Domains used for Uniprot and UniRef
my @uni_domains = ('ftp://ftp.ebi.ac.uk', 'ftp://ftp.expasy.org', 'ftp://ftp.uniprot.org');
###########################################


=head2 check_url

 Title   : check_url
 Usage   : check_url(link => the link to check);
 Function: Checks whether the link is working by sending a HEAD request with curl. Supports both HTTP and FTP links.
           Tries to probe the URL 5 times before declaring failure.
 Returns : 1 if link is good, 0 otherwise

=cut
sub check_url {
    my $url = shift;
    my $attempts = 0;
    while ($attempts < 5) {
        # -f means die on 404 and similar errors
        # -s means do not print the progress bar
        # --head specified type of request (HEAD)
        # --location tells curl to follow redirects, if any
        system "curl -f -s --head --location $url >/dev/null";
        if (!$?) {  # All good, link is alive
            return 1
        }
        $attempts = 1;
        sleep(1);  # trying again *right away* is not very wise, let's wait at least a little
    }
    # If we're here, all attempts ended in failure, link is dead
    return 0
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
      }else{
     
        #This command will force the download with different methods
        if (ftp_fetch_file($fileAdd,$folderName,$user,$pass) == 1){
            print "..completed!\n";
            $retVal = 1;
        #If it does not work and we are searching Uniprot, then access different servers
        }elsif (match_array_with_string($fileAdd, \@uni_domains) == 1){
          $new_link  = download_from_diff_sources($fileAdd,$folderName,$user,$pass);
          $retVal = 2;
        }
      }
  }
  if ($retVal == 0){
    die "Unable to download file $fileAdd. Restart again later. Annocript will close...\n";
  }
  return $retVal, $new_link;
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
			 #Download
       #print "\n Downloading File: ".$link."\n";
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
		#Uncompression
    #print "\n Uncompress File: ".$dbDataFolder."/".$compressed."-->".$dbDataFolder."/".$uncompressed."\n";
    my_extract_any_file($dbDataFolder."/".$compressed,$dbDataFolder);
 }else{print "File ".$dbDataFolder."/".$uncompressed." already present...\n";}#Exists the file 
}


=head2 download_from_diff_sources 

 Title   : download_from_diff_sources
 Usage   : download_from_diff_sources( -link => file to download from different sources
                               );

 Function: This function only works when the database is about Uniprot. It searches the files on EBI and Expasy if
 Uniprot is not working and viceversa
 
 Returns : 1 if succeed in changing the link, -1 otherwise. And new_link is the new link to change in the config hash

=cut
sub download_from_diff_sources{
  
   my $link = shift;
   my $folderName = shift;
   my $user = shift;
   my $pass = shift;
   
   my $new_link = '';
   my $retVal = 0;
   
   my @domains = @uni_domains;
   
   #Remove the link already tried 
   foreach my $domain (@domains){
     if ($link =~ /$domain/){
        @domains = grep { $_ ne $domain } @domains;#Remove the corresponding element from the array
      }
    }
    
   #If the /pub is there, it means that the link is with Uniprot.org or EBI. 
   #Thus add it for expasy and not for the other
   foreach my $domain (@domains){
     if ( ($domain =~ /\/pub\//) and ($retVal == 0)){#($link =~ /\/pub\//) and 
       my @link_pieces = split("/", $link);#separate the link pieces
       my @dom2_pieces = split("/",$domain); #separate the pieces in the other domain
       splice(@link_pieces,0,3,@dom2_pieces);#substitute the dom2 at the dom1
       $new_link = join("/",@link_pieces);#new link is composed by joining!
       
       #Now access it to download!
       if(ftp_fetch_file($new_link,$folderName,$user,$pass) == 1) {
        $retVal = 1;
        print "Downloaded $new_link because $link was not reacheble...\n";
       }
     }
     #If the dom1 is with /pub and the second is not. Then substitute differently differently
     if ( ($domain !~ /\/pub\//) and ($retVal == 0)){#($link =~ /\/pub\//) and 
       my @link_pieces = split("/", $link);#separate the link pieces
       my @dom2_pieces = split("/",$domain); #separate the pieces in the other domain
       splice(@link_pieces,0,4,@dom2_pieces);#substitute the dom2 at the dom1
       $new_link = join("/",@link_pieces);#new link is composed by joining!
       
       #Now access it to download!
       if(ftp_fetch_file($new_link,$folderName,$user,$pass) == 1) {
        $retVal = 1;
        print "Downloaded $new_link because $link was not reacheble...\n";
       }
     }
    }
    
   #if ($retVal == 0){
      #delete_file($folderName."/".$fileToDl);
      #die "ERROR [$?]: Annocript can't download $link. Servers are all out! Please try again later.\n"; 
      #$retVal = -1;
   #}
    
  return $new_link;
}



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
   $File::Fetch::WARN = 0;
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
    my $new_link = '';
    
    my $ebi_dom = 'ftp://ftp.ebi.ac.uk';
    my $exp_dom = 'ftp://ftp.expasy.org';
    
    my $req = HTTP::Request->new(GET => $link);
    my $UA = new LWP::UserAgent;
    my $res = $UA->request($req);
       
    if($res->is_success) {
      print $res->as_string();
      #$retVal = $res->header('Content-Length');
      $retVal = 1;
    } else {
          #If does not work and is uniprot, we can try with EBI
          if ( (split(/ /,$res->status_line))[0] eq '404'){
            print "WARNING:(404 File doesn't exists at $link). If you are sure of this link leave Annocript running.\n";
            if ($link =~ /\/pub\// and $link =~ /uniprot.org/){
              my @link_pieces = split("/", $link);
              my @ebi_pieces = split("/",$ebi_dom);
              splice(@link_pieces,0,3,@ebi_pieces);
              $new_link = join("/",@link_pieces);
              print "Trying to access with EBI link $new_link\n";
              print "Request in corso...\n";
              if(head($new_link)){
                $retVal = 2;
              }else{
                #If EBI also does not work, Try with expasy
                my @link_pieces = split("/", $link);
                my @exp_pieces = split("/",$ebi_dom);
                splice(@link_pieces,0,4,@exp_pieces);
                $new_link = join("/",@link_pieces);
                print "Trying to download with EBI link $new_link\n";
                print "Request in corso...\n";
                if(head($new_link)){ 
                  $retVal = 2;
                }else{
                  $retVal = -1;
                }
               }
            }else{
              #Itis Expasy link. Try with EBI.
              my @link_pieces = split("/", $link);
              my @ebi_pieces = split("/",$ebi_dom);
              splice(@link_pieces,0,3,@ebi_pieces);
              $new_link = join("/",@link_pieces);
              print "Trying to download with EBI link $new_link\n";
              if(head($new_link)){
                $retVal = 2;
              }else{
                $retVal = -1;
              }
            }
          }
      }      
  return $retVal, $new_link;
}


=head2 testFTP

 Title   : testFTP
 Usage   : testFTP( - link => the link to check
                    );

 Function: It checks if a link to an FTP external resource is working. If it is it returns 1. If it is not then it will check
          what kind of error is coming out and alert the user with a message.

 Returns : 1 if link is good, 0 otherwise

=cut
sub testFTP {
    my $ftpLink = shift;
    my $user = shift;
    my $pass = shift;
    
    my ( $host, $filepath ) = $ftpLink =~ m|ftp://([^/]+)(.*)$|;
    
    my $status= 0;
    my $numTent = 0;
    
    $user = 'anonymous' unless defined $user;
    $pass =  'nobody@nowhere.com' unless defined $pass;
    
    while ( ($numTent < 5) and ($status != 1)){
        if (my $ftp = Net::FTP->new( $host )) {
          $ftp->login(  $user, $pass );
          my $size = $ftp->size( $filepath );
          $ftp->quit;
          $status = defined $size ? 1 : -1;
        }
        else {
            $status = -1;
        }
        $numTent++;
    }
    
    
  return $status;
}


=head2 check_FTP_diff_sources 

 Title   : check_FTP_diff_sources
 Usage   : check_FTP_diff_sources( -link => file to download from different sources
                               );

 Function: This function tests the FTP of different domains for Uniprot databases and if some of them is working it gives
  in output a new link to substitute the actual
  
 Returns : 1 if succeed in changing the link, -1 otherwise. And new_link is the new link to change in the config hash

=cut
sub check_FTP_diff_sources{
  
   my $link = shift;
   my $user = shift;
   my $pass = shift;
   
   my $new_link = '';
   my $retVal = 0;
   
   my @domains = @uni_domains;
   
   #Remove the link already tried 
   foreach my $domain (@domains){
     if ($link =~ /$domain/){
        @domains = grep { $_ ne $domain } @domains;#Remove the corresponding element from the array
      }
    }
    
   #If the /pub is there, it means that the link is with Uniprot.org or EBI. 
   #Thus add it for expasy and not for the other
   foreach my $domain (@domains){
     if ( ($domain =~ /\/pub\//) and ($retVal == 0)){#($link =~ /\/pub\//) and 
       my @link_pieces = split("/", $link);#separate the link pieces
       my @dom2_pieces = split("/",$domain); #separate the pieces in the other domain
       splice(@link_pieces,0,3,@dom2_pieces);#substitute the dom2 at the dom1
       $new_link = join("/",@link_pieces);#new link is composed by joining!
       
       #Now access it to download!
       if(testFTP($new_link,$user,$pass) == 1) {
        $retVal = 1;
        print "Annocript will use $new_link because $link was not reachable...\n";
       }
     }
     #If the dom1 is with /pub and the second is not. Then substitute differently differently
     if ( ($domain !~ /\/pub\//) and ($retVal == 0)){#($link =~ /\/pub\//) and 
       my @link_pieces = split("/", $link);#separate the link pieces
       my @dom2_pieces = split("/",$domain); #separate the pieces in the other domain
       splice(@link_pieces,0,4,@dom2_pieces);#substitute the dom2 at the dom1
       $new_link = join("/",@link_pieces);#new link is composed by joining!
       
       #Now access it to download!
       if(testFTP($new_link,$user,$pass) == 1) {
        $retVal = 1;
        print "Annocript will use $new_link because $link was not reachable...\n";
       }
     }
    }
    
   #if ($retVal == 0){
      #delete_file($folderName."/".$fileToDl);
      #die "ERROR [$?]: Annocript can't download $link. Servers are all out! Please try again later.\n"; 
      #$retVal = -1;
   #}
    
  return $new_link;
}






=head2 paths_hash_to_html

 Title   : paths_hash_to_html
 Usage   : paths_hash_to_html( config_file => the config hash );

 Function: Reads the paths hash obtained from the outlist file and writes a file in HTML
					with for each PAGE and for each SAMPLE and for each READFILE the description linking to 
					the given path into the webserver.
					The variable $host is used to change the host name 
  
 Returns : a string that is actually an HTML page
 
=cut
sub paths_hash_to_html{
	my $cfg_hash = shift;
	my $page = shift;
	my $paths_hash = shift;
	my $host = shift;
	my $workingFolder = shift;
		
	my $cgi = CGI->new; #Instantiate a CGI class
	
	my $string = "<center>";

		#
		while (my ($sample_name, $value) = each %{ $paths_hash->{$page} } ) {
			#Print the sample name if exists
			if ($sample_name ne '-'){
				$string = $string."Sample: $sample_name<br>";
			}
			while (my ($readf_name, $value2) = each %{ $paths_hash->{$page}->{$sample_name}} ) {
				while (my ($desc, $value3) = each %{ $paths_hash->{$page}->{$sample_name}->{$readf_name}} ) {
					#if ($desc ne '-'){print FILE $desc."\n"}
					#Extract the path and make the link
					my $from_wf = $paths_hash->{$page}->{$sample_name}->{$readf_name}->{$desc}->{'path'};
					if (defined $from_wf) {
						$from_wf =~ s/$workingFolder//;
						my $link = $host.$from_wf;
						#Concatenate description with sample and read file name
						if ($sample_name ne '-'){$desc = "$desc $sample_name ";}
						if ($readf_name ne '-'){$desc = "$desc $readf_name ";}
						#Concatenate the link
						$string = $string.$cgi->a({href=>$link},$desc)."<br>";						#code
					}
					

				}
			}
		}

 $string .= "</center>";
 #print_and_log( "String is (inside) $string ..\n",$log_file);
 return $string;
}

=head2 file_2_html

 Title   : file_2_html
 Usage   : file_2_html(   );

 Function: Reads a file line by line and transforms it in a TAGged HTML text
					It is supposed that each line is a distinct paragraph
 Returns : a string with the HTML code
 
=cut
sub file_2_html{
	my $inFile = shift;

 #Create
 #print_and_log( "Opening $inFile ..\n",$log_file);
 open(INFILE, "<$inFile") or warn "WARNING [$!]: can't open file: $inFile. Check permissions.\n"; 
 my $string = "<center>";
 while (my $line = <INFILE>){
		chop($line);
		#Use <p> to separate paragraphs
		$string = $string."<p>".$line."</p>\n";
		#print_and_log( "Line read $line ..\n",$log_file);
 }
 close(INFILE);
 $string .= "</center>";
 #print_and_log( "String is (inside) $string ..\n",$log_file);
 return $string;
}


=head2 table_list_2_html

 Title   : table_list_2_html
 Usage   : table_list_2_html(   );

 Function: Transforms a list of tables in HTML. The input is a file with a list of
						TableName TablePath
						And the tables are text file tab separated. The function table_2_html is
						iteratively called.
					
 Returns : a string with the HTML code
 
=cut
sub table_list_2_html{
	my $tables_list_file = shift;
	my $remove_temp = shift;

	my $html_string = "";

	
	if (-e $tables_list_file) {
		#print "Opening $tables_list_file ..\n";
		open(INFILE, "<$tables_list_file") or print "WARNING [$!]: can't open file: $tables_list_file. Check permissions.\n"; 
		
		#Transform all the tables in html
		my @tables_paths = list_to_array($tables_list_file,'NO_NEW_LINE');
		foreach my $table_path (@tables_paths){
			my @fields = split("\t",$table_path);
				if ( -e $fields[1]){
				$html_string .= "<p>".$fields[0];
				$html_string .= table_2_html($fields[1]);
				$html_string .= "</p>\n";
				if ( $remove_temp eq 'YES' ){
					delete_file($fields[1]);
				}
		 }
		}
		close(INFILE);		
	}#else{print "WARNING [$!]: file $tables_list_file does not exist.\n";}

	return $html_string;
}


=head2 table_2_html

 Title   : table_2_html
 Usage   : table_2_html(   );

 Function: Transforms a file with a table (columns tab separated)
					in a TAGged HTML text
					N.B. First line contains the header of the table! columns are TAB separated!!
					
 Returns : a string with the HTML code
 
=cut
sub table_2_html{
	my $inFile = shift;

	#print_and_log( "Opening $inFile ..\n",$log_file);
	open(INFILE, "<$inFile") or die "ERROR [$!]: can't open file: $inFile. Check permissions.\n"; 
	my $string = '<center><table border="1" style="width:50%">';

	my $num_row = 0;
	while (my $line = <INFILE>){
		chop($line);
		my @fields = split("\t",$line);
		$string .= "\t<tr>\n";
		#Header
		if ($num_row == 0){
			foreach my $field (@fields){
					$string .= "\t\t<th>$field</th>\n";
			}
		}else{
			foreach my $field (@fields){
					$string .= "\t\t<td>$field</td>\n";
			}			
		}
		$string .= "\t</tr>\n";
		$num_row++;
		#print_and_log( "Line read $line ..\n",$log_file);
	}
	close(INFILE);
	$string .= "</table></center>\n";
	#print_and_log( "String is (inside) $string ..\n",$log_file);
	return $string;
}


=head2 var_stats_2_html

 Title   : var_stats_2_html
 Usage   : var_stats_2_html(  config_file => the config hash
								);

 Function: Prints an HTML page using the statistics generated with R
					for the result of VarGenius.
 
 Returns :
 
=cut
sub var_stats_2_html{
	my $var_stats_f = shift;
	
	my $datasep = shift;
	
	my $totreads_block_ind = 1;
	my $freqs_block_ind = 2;
	my $segregation_block_ind = 3;
	my $dbsnpexist_block_ind = 4;
	my $zygosity_block_ind = 5;
	my $zyglegenda_block_ind = 6;
	
	#print_and_log( "Opening $inFile ..\n",$log_file);
	 open(INFILE, "<$var_stats_f") or die "ERROR [$!]: can't open file: $var_stats_f. Check permissions.\n"; 
	 my $string = '<center>';
	 my $start = 0;
	 my $tab_line = 0;
	 while (my $line = <INFILE>){
			chop($line);
			#A block of information starts with a set of # and a string
			if ($line =~ /#+\w+/){
				$start++;
				$line =~ s/#//g;
				$string .= "<br>".$line."<br>";
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind) {
					$tab_line = 0;
					$string .= '<table border="1" style="width:50%">';
				}
				#$line = <INFILE>; 
			}#and ends with only a set of #
			elsif ($line =~ /#+/) {
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind){
					$string .= '</table>';
				}
				$string .= "<br>";
			}else{
				#This instructions are to read and build tables
				if ( $start == $totreads_block_ind or $start == $zygosity_block_ind or $start == $segregation_block_ind){
					$string .= "<tr>";
					my $tag_s;
					my $tag_e;
					my @parts = split("\t",$line);		
					if ($tab_line == 0){
						$tag_s = "<th>";
						$tag_e = 	"</th>";
						$tab_line++;
					}else{
						$tag_s = "<td>";
						$tag_e = 	"</td>";					
					}
					for (my $i = 0; $i <scalar(@parts) ; $i++){
							$string .= $tag_s.$parts[$i].$tag_e;
					}
					$string .= "</tr>";	
				}

				#Second table
				if ( $start == $freqs_block_ind){
					#if ($line =~ /^Freq/){
						#$string .= "<br><br>".$line."<br>";
						##$line = <INFILE>; 
					#}els
					if ($line =~ /^Chrom/){
						$string .= "<br>".$line." : ";
						#$line = <INFILE>; 
					}elsif ($line ne ''){
						$string .= $line.$datasep;
					}
				}				
	
				#Legenda
				if ( $start == $zyglegenda_block_ind or $start == $dbsnpexist_block_ind){
						$string = $string."<h6>".$line."</h6>";
				}					
			}
				

			
			#print_and_log( "Line read $line ..\n",$log_file);
	 }
	 close(INFILE);
	 $string .= "</center>";
	 return $string;
}


=head2 html_page

 Title : html_page
 Usage : html_page(  imgDir -> #Directory where the images are stored	
					 htmlFile -> Path to the html file to create
					 ($descrText) -> array pointer descriptive texts of the links to the other pages
					  $descrTextIndex -> Index of the descriptive text of this page
					 ($links) -> pointer to array with the links at the other pages
					 ($imgsKeywords)-> pointer to array of keywords to use to search the images inside the folder
		);

 Function: creates an html page using a directory where are searched the images. Prints the header using $descrTextIndex
			from ($descrText). Puts the logo from the path. A loop searches for all the images in $imgDir
			with title comprising a word inside ($imgsKeywords). Print a series of links to the other pages using a loop
			searching in ($links) whose names are in descrText. Finally prints information at the bottom of the page.

 Returns : nothing

=cut 
sub html_page{
 my $cfg_hash = shift;
 my $imgPath = shift;#Directory where the images are stored	
 my ($descrText) = shift;#descriptive texts of the links to the other pages
 my $descrTextIndex = shift;#Index of the descriptive text of this page
 my ($links) = shift;#array with the links at the other pages
 my $imgsKeywords = shift;#keywords to use to search the images inside the folder
 my $infoText = shift;#A string to print as informations in the BODY of the page
 my $html_folder = shift;
 
 my $cgi = CGI->new; #Instantiate a CGI class
 
 #Open a new html file
 open (OUT,">$html_folder/".$$links[$descrTextIndex]) || die "ERROR [$!]: can't open file ".$$links[$descrTextIndex]."\n";
  
 print OUT $cgi->start_html("VarGenius ".$cfg_hash->{'VarGeniusVer'}."- ".$$descrText[$descrTextIndex]), # Begin HTML page
       #$cgi->header, # First we create a header
       $cgi->center; 
 if ( $imgPath ne 'none'){
	print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$cfg_hash->{'VarGeniusLogo'}."' ALT='image' ALIGN=TOP> <br>";
       #$cgi->h1('<i>Annocript '.$configHash->{'AnnocriptVer'}.'</i>'), # create the Tag <h1>TEXT</h1>
 }
 print OUT $cgi->h3($$descrText[$descrTextIndex]."\n");
 #BODY
 print OUT $cgi->p($infoText);
 
 #Open images folder if it exists
 if ( $imgPath ne 'none'){
	 opendir IMG, $imgPath or die "ERROR [$!]: cannot open dir $imgPath.\n";
	 my @images= readdir IMG;
	 my @imgsSorted = sort {$a cmp $b} @images;#order alfabetically the images names
	closedir IMG;
	 #Search and write HTML code to show images if it exitsts
	 #print OUT "</center>";
	 if(defined($imgsKeywords)){
		 foreach my $img (@imgsSorted){
			#if (scalar(grep {/\b$img\b/} @$imgsKeywords) == 1){
			if ($img =~ /$imgsKeywords/){
				print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$img."' ALT='image' ALIGN=TOP> <br>";
			} 
		 }
	 }
 }
 

	#FOOTER
 #print the website map
 for ( my $link=0; $link<scalar(@$links);$link++ ){
	print OUT $cgi->a({href=>$$links[$link]},$$descrText[$link]);
	print OUT " | ";
 }

 #Print information
 print OUT $cgi->h5('VarGenius '.$cfg_hash->{'VarGeniusVer'}.' - Copyright of Tigem'),
      $cgi->p(scalar localtime); # create <p>Text</p>
 print OUT "</center>";
 print OUT $cgi->end_html; # end of HTML
 close(OUT);
}


1;
