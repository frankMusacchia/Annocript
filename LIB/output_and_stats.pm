#############################output_and_stats.pm#############################################
# This module is to extract information from the output of Annocript and create some plots  #
# and a web page with statistics.                                                           #
#############################################################################################

#Annocript - A complete tool for transcriptomes annotation
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

package LIB::output_and_stats;

use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( lncRNA_coding_fasta_file_generator lncRNA_fasta_file_generator 
extractAnnotated extractStatistics );

#GENERAL UTILITY MODULES
use strict;
use warnings;
use Data::Dumper;#To print the dumper of the hashes
use Cwd;#get pathname of current working directory
use CGI; #load cgi module to produce HTML output
use List::MoreUtils qw{ any all }; #Used to search in arrays
use Time::HiRes qw( time );
use IO::Handle;#To immediately print with autoflush 

#FILE MANAGEMENT MODULES
use File::Copy; #Used to move files


#Using a library for standard utilities								
use LIB::std_lib qw(print_array);	

#Using a library to manage files
use LIB::files_management qw(  extract_name checkOutputsPresence check_presence);


use LIB::annocript_management qw();

use LIB::html_management qw(dl_and_extract download_file);
									
									
##DECLARATION OF GLOBAL INDEXES OF THE FIELDS IN THE ANNOCRIPT OUTPUT
my $transcrId ; my $transcrLength ;
#my $count;
my $Name1; my $Length1;
my $Score1; my $HitLength1;
my $QCoverage1; my $HCoverage1; 
my $Strand1;
my $Frame1;
my $Descr1;
my $closerOS1;
my $EnzId; my $EnzDescr; 
my $PwLev1;my $PwLev2;my $PwLev3;

my $Name2 ; my $Length2; 
my $Score2; my $HitLength2;
my $QCoverage2 ; my $HCoverage2;
my $Strand2;
my $Frame2;
my $Descr2;
my $closerOS2;

my $BPId; my $BPDesc;
my $MFId; my $MFDesc;
my $CCId; my $CCDesc;

my $domBPId; my $domBPDesc;
my $domMFId; my $domMFDesc;
my $domCCId; my $domCCDesc;

my $CDName; my $CDStartEnd;
my $CDScore ; my $CDDescr;

my $riboName; my $riboScore;
my $riboDesc;
my $longORFLength; my $longORFS; my $longORFF;

my $probNC; my $NC4Annocript; my $seq;


######GLOBAL VARIABLES USED
my $configHash;#A variable for the configuration hash
my $dsn; #Used for database configuration
my $numWords = 10;#Number of words to use for write descriptions
my $maxLengthDescs;#Maximum Number of chars for the descriptions in files and plots
my $separator;#The separator used
my $database;#Database used
my $cdName4Expression;#The Conserved Domains database to use for plots
my $outFolder = '';#Path of the folder with the outputs
my $logFolder = '';#Path of the folder with the log files
my $statsFolder = '';
my $headHash; #An hash which contains the reference to all the header fields

################################################################################################

=head2 settingsForUniref

 Title : settingsForUniref
 Usage : settingsForUniref(  );

 Function: this function sets the indexes of the array  containing Annocript output for the case of using the
            uniref

 Returns : nothing

=cut
sub setFieldNames {
  my $outFileName = shift;
  my $dbInUse = shift;
  
  #Open the output file and read the header
  open(OUT_ANN,"<$outFileName") or die "Cannot open the file $outFileName: ?\n";
  my $row = <OUT_ANN>;
  chomp($row);
  #print "Just read the first line of $outFileName: ".$row;
  #Put the headers name in an hash and associate them the field number
  my @fields = split(/\t/, $row);
  #print_array(\@fields);#DEBUGCODE
  #print "Number of fields :".scalar(@fields)."\n";#DEBUGCODE
  my $numField = 0;
  #Here we fill the hash with the field's names and assign a number
  foreach my $field (@fields){
	  $headHash->{$field} = $numField;
	  #print "Assigning $numField to the field '$field': ".$headHash->{$field}."\n";
	  $numField++;
	  
  }	
  close(OUT_ANN);
    
  ##INDEXES OF THE FIELDS IN THE ANNOCRIPT OUTPUT
  ## IF THE FIELDS CHANGE THEN YOU MUST CHANGE THE NAMES ACCORDINGLY
  ## SEE ALSO HEADERS IN GFF3_MANAGER.pm
  $transcrId = 'TranscriptName'; $transcrLength = 'TransLength';

  $Name1 = 'HSPNameSP'; $Length1 = 'HSPLengthSP'; 
  $Score1 = 'HSPEvalueSP'; $HitLength1 = 'HITLengthSP';
  $QCoverage1 = 'QCoverageSP';  $HCoverage1 = 'HCoverageSP'; 
  $Strand1 = 'StrandSP';
  $Frame1 = 'FrameSP';
  $Descr1 = 'DescriptionSP'; 
  $closerOS1 = 'OSNameSP';  
  $EnzId = 'EnzymeIds';  $EnzDescr = 'EnzymeDescs'; 
  $PwLev1 = 'PwLev1'; $PwLev2 = 'PwLev2'; $PwLev3 = 'PwLev3';
  
  if ($dbInUse eq "uniprotkb"){ 
		$Name2 = 'HSPNameTR'; $Length2 = 'HSPLengthTR' ; 
		$Score2 = 'HSPEvalueTR'; $HitLength2 = 'HITLengthTR';
		$QCoverage2 = 'QCoverageTR';  $HCoverage2 = 'HCoverageTR'; 
		$Strand2 = 'StrandTR';
		$Frame2 = 'FrameTR';
		$Descr2 = 'DescriptionTR'; 
		$closerOS2 = 'OSNameTR';
  }
  if ($dbInUse eq "uniref"){ 
		$Name2 = 'HSPNameUf'; $Length2 = 'HSPLengthUf' ; 
		$Score2 = 'HSPEvalueUf'; $HitLength2 = 'HITLengthUf';
		$QCoverage2 = 'QCoverageUf';  $HCoverage2 = 'HCoverageUf'; 
		$Strand2 = 'StrandUf';
		$Frame2 = 'FrameUf';
		$Descr2 = 'DescriptionUf'; 
		$closerOS2 = 'Taxonomy';
	}   
 
  $CDName = 'CDName'; 
  $CDStartEnd = 'CDStartEnd';
  $CDScore = 'CDEvalue';  $CDDescr = 'CDDesc';

  $BPId = 'BPId';  $BPDesc = 'BPDesc';
  $MFId = 'MFId';  $MFDesc = 'MFDesc';
  $CCId = 'CCId';  $CCDesc = 'CCDesc';

  $domBPId = 'domBPId';  $domBPDesc = 'domBPDesc';
  $domMFId = 'domMFId';  $domMFDesc = 'domMFDesc';
  $domCCId = 'domCCId';  $domCCDesc = 'domCCDesc';

  $riboName = 'OtherNCName'; $riboScore = 'OtherNCEvalue'; $riboDesc = 'OtherNCDesc'; 
  $longORFLength = 'LongOrfLength';
  $longORFS = 'LongOrfStrand';  $longORFF = 'LongOrfFrame';

  $probNC = 'ProbToBeNonCoding';  $NC4Annocript = 'lncRNA4Annocript'; $seq = 'Sequence';

}


=head2 lncRNA_fasta_file_generator

 Title : lncRNA_fasta_file_generator
 Usage : lncRNA_fasta_file_generator( -annocriptOut => 'the output of Annocript',
         -outName => 'the name of the output file',
         -outFolder => 'the output folder'
     );

 Function: this function extract from the Annocript's output the non coding sequences by using the corresponding field

 Returns : nothing

=cut
sub lncRNA_fasta_file_generator {
 
 my $annocriptOut = shift;
 my $outName = shift;
 my $outFolder = shift;
 
 my $NCCount = 0;
 $annocriptOut = $outFolder.'/'.$annocriptOut; 
 #print "The path to ann out: $annocriptOut\n";
 print "Creating a file with the long non-coding RNA sequences selected by Annocript....";
 if (open (ANNOUT_FILE,"<$annocriptOut")){
  open (NONCOD_FILE,">$outName");
  

  #Creation of a file with the non coding sequences
  my $row = <ANNOUT_FILE>; #First line to not to read
  while( $row = <ANNOUT_FILE> ) {
    my @fields = split ("\t",$row);
    if ( $fields[$headHash->{$NC4Annocript}] eq '1'){
      #..then it is non coding and write it in a file
      print NONCOD_FILE ">".$fields[$headHash->{$transcrId}]."|".$fields[$headHash->{$longORFLength}]."|".$fields[$headHash->{$probNC}]."\n";
      my $lncRNAseq = "";
      while ($fields[$headHash->{$seq}] =~ /(.{1,60})/g) {
				$lncRNAseq= $lncRNAseq.$1."\n";
			}
      print NONCOD_FILE $lncRNAseq;
      #print NONCOD_FILE join("\t",$fields[$headHash->{$transcrId}],$fields[$headHash->{$longORFLength}],$fields[$headHash->{$probNC}],$fields[$headHash->{$seq}],"\n");
      $NCCount++;
    } 
  }
  close (ANNOUT_FILE);
  close (NONCOD_FILE);
  move($outName,$outFolder) or print "ERROR [$!]: unable to move $outName in $outFolder. Please check permissions.\n";#Moving the output of lncRNA predictor to the output folder
 }else {print "$annocriptOut is missing. Something has been wrong!\n";}
 print "..DONE!\nThe output has been wirtten in $outFolder/$outName \n";
 return $NCCount;
}
 

=head2 lncRNA_coding_fasta_file_generator

 Title : lncRNA_coding_fasta_file_generator
 Usage : lncRNA_coding_fasta_file_generator( -annocriptOut => 'the output of Annocript',
         -outName => 'the name of the output file',
         -outFolder => 'the output folder'
     );

 Function: this function extract from the Annocript's output the non-coding and coding sequences and creates
      two corresponding fasta files.

 Returns : nothing

=cut
sub lncRNA_coding_fasta_file_generator {
 
 my $annocriptOut = shift;
 my $lncRNAFileName = shift;
 my $codingFileName = shift;
 my $outFolder = shift;
 my $NCORFLength = shift;
 
 
 my $NCCount = 0;
 my $itsCoding = 0;
 $annocriptOut = $outFolder.'/'.$annocriptOut; 
 print "The path to ann out: $annocriptOut\n";
 if (open (ANNOUT_FILE,"<$annocriptOut")){
  open (NONCOD_FILE,">$lncRNAFileName") or die "Cannot open $lncRNAFileName $?";
  open (CODING_FILE,">$codingFileName") or die "Cannot open $codingFileName $^";
  

  #Creation of a file with the non coding sequences
  my $row = <ANNOUT_FILE>; #First line not to read
  while( $row = <ANNOUT_FILE> ) {
    my @fields = split ("\t",$row);
    if ( $fields[$headHash->{$NC4Annocript}] eq '1'){
      #..then it is non coding and write it in a file
      print NONCOD_FILE ">".$fields[$headHash->{$transcrId}]."|".$fields[$headHash->{$longORFLength}]."|".$fields[$headHash->{$probNC}]."\n";
      my $lncRNAseq = "";
      while ($fields[$headHash->{$seq}] =~ /(.{1,60})/g) {
		$lncRNAseq= $lncRNAseq.$1."\n";
	  }
      print NONCOD_FILE $lncRNAseq;
      #print NONCOD_FILE join("\t",$fields[$headHash->{$transcrId}],$fields[$headHash->{$longORFLength}],$fields[$headHash->{$probNC}],$fields[$headHash->{$seq}],"\n");
      $NCCount++;
    }#otherwise the coding sequence will be written in a different file
      elsif ( ($fields[$headHash->{$Name1}]  ne '-') or ( $fields[$headHash->{$Name2}] ne '-') 
					or  ($fields[$headHash->{$CDName}] ne '-')  ){
				$itsCoding = 1;
      }elsif ( $fields[$headHash->{$longORFLength}] ne '-' ){#Avoid to meet a not analyzed value
					if ( $fields[$headHash->{$longORFLength}] > $NCORFLength){
						$itsCoding = 1;
					}
				} 
		if ( $itsCoding == 1){
			#..then it is coding and write it in a file
			print CODING_FILE ">".$fields[$headHash->{$transcrId}]."|".$fields[$headHash->{$longORFLength}]."|".$fields[$headHash->{$probNC}]."\n";
			my $codingSeq = "";
			#print "$seq sta nell'indice:".$headHash->{$seq}."\n";
			while ($fields[$headHash->{$seq}] =~ /(.{1,60})/g) {
				$codingSeq= $codingSeq.$1."\n";
			}
			print CODING_FILE $codingSeq;
			$itsCoding = 0;		
		}
  }
  close (ANNOUT_FILE);
  close (NONCOD_FILE);
  close (CODING_FILE);
  move($lncRNAFileName,$outFolder) or print "ERROR [$!]: unable to move $lncRNAFileName in $outFolder. Please check permissions.\n";#Moving the lncRNA file in the output folder
  move($codingFileName,$outFolder) or print "ERROR [$!]: unable to move $codingFileName in $outFolder. Please check permissions.\n";#Moving the coding seqs file
 }else {print "$annocriptOut is missing. Something has been wrong!\n";}
 print "..DONE!\nFiles with coding sequences and lncRNA have been put in $outFolder\n";
 return $NCCount;
}

=head2 extractGOstats

 Title : extractGOstats
 Usage : extractGOstats( -dbh => 'Reference to the database accessed',
                         -outAnnArray => 'the output of Annocript in a matrix',
                         -iRow => 'the current row analyzed',
                         -id => 'Number of the field in the table of the go ID'
                         -hash => 'The GO hash that we are building with the statistics'
                         
     );

 Function: this function builds an hash counting the occurrence of a given GO term and its expression 
 level using the counts

 Returns : nothing but the $hash inserted will be modified

=cut
sub extractGOstats {
  #Input
   my $dbh = shift;#Reference to the database accessed
   my $outAnnArray = shift;#The matrix with all the table of output of Annocript
   my $iRow =shift; #the current row analyzed
   my $id = shift;#Number of the field in the table of the ID
   my ($hash) = shift;#The GO hash that we are building with the statistics
   my $countsPresent = shift;#Says if the counts are present for the enrichments
   my $association = shift;
   
   #Local variables
   my @insertedEls = ();#Elements in the row that have been already inserted
   my $splitInd = 0;
   my @splitted;
   
   #For domains separates the GO terms separated by semicolon, too 
   if (defined $association){
	   if($association eq 'domains'){
		  my @partSplitted = split(/\Q$separator\E/,${$outAnnArray}[$iRow][$id]);#Array to split ids
		  foreach my $singleSplit (@partSplitted){
		    push(@splitted,split(/;/,$singleSplit));	  
		  }
		}
   }else{
	  @splitted = split(/\Q$separator\E/,${$outAnnArray}[$iRow][$id]);#Array to split ids 
   }
   
   #The loop goes for each of the elements separated
   foreach my $goTerm (@splitted){
     #and if the element hasn't been already evaluated
     if ( scalar(grep {/\b$goTerm\b/} @insertedEls) == 0 ){
        push(@insertedEls,$goTerm);#Insert the go term in the array of already taken ones
        
        #push it in the already examined list and access db to fetch its informations
        my @res = ();#Result of the query
        my $select_handle;
        #Extracts information from UNIPROT only if $database is defined. (blastx i.e.)
        if (defined $database){
          #fetch informations from the database. Using the database will not bring to possible incorrect informations coming from
          #the output files. 
          $select_handle = $dbh->prepare ("SELECT division, definition FROM goTable WHERE go_id='".$goTerm."'");
          die "ERROR [$?]: can't prepare queries; aborting: ?" unless defined $select_handle;
          $select_handle->execute() or die $select_handle->errstr;
          @res = $select_handle->fetchrow_array;
        }else{print "ERROR in output_and_stats: $database was not connected in a point where Annocript don't expect\n";}
        
        #If results are coming out from the database...
        if ( scalar(@res) > 0 ){
         #..we take the description from the result and separate it to get only $maxLengthDescs chars.
         #this is done to not have a description too long (it will be printed on an image with R).
         my $desc = $res[1];
        
         my $descString = substr $desc, 0 ,$maxLengthDescs;#Reducing the dimensions of description
         
         $$hash->{$goTerm}->{'count'}++;
         $$hash->{$goTerm}->{'desc'} = $descString;
         
        }else{print "ERROR in output_and_stats: for the GO Term: $goTerm has not been found a corresponding term in the database.\n";}
     }else{
			 if (not defined $association){
				print "ERROR in output_and_stats: two equal GO terms have been found on line $iRow of the output file. ".
							"Some error is present in the creation in gff_manager module. Please contact".
              " the creator of Annocript.\n";}
			 }
   }

}




=head2 extractStatistics

 Title : extractStatistics
 Usage : extractStatistics( -configHash => 'the hash with all the variables for the configuration',
                          -trWords => 'an array with the words in TREMBL to be considered non-coding',
     );

 Function: this function extract from the Annocript's output a set of statistics:
    - mean length of sequences
    - minimum, maximum sequences length
    - Mean Percentage of Adenine, Guanine, Thymine, Cytosine, N
    - Number of annotated sequences
    - Mean Percentage of GC content
    
    It also create files for:
     - gc content :sequence -> gccontent
     - biological process: GOID -> percentage of sequences referred to it
     - molecular function: GOID -> percentage of sequences referred to it
     - cellular component: GOID -> percentage of sequences referred to it 
     - domains: domain -> number of sequences referred to it
     
  it uses the uniprot DB to fetch informations.
    
 Returns : a file with statistics - stats_name of the query file.fasta

=cut
sub extractStatistics{
 
 $configHash = shift;
 my @trWords = shift;#Words to match to be a non coding
 
 STDOUT->autoflush(1);#This makes STDOUT hot in the sense that everything will be print immediately
 
 #Import variables from configuration file 
 my $annocriptOutFile = $configHash->{'outFiltered'};#Output file of Annocript
 #my $annocriptOutFile = $configHash->{'outFileName'};#Output file of Annocript
 my $querySeqs = $configHash->{'fastaSeqs'};#fasta file with sequences
 my $orfFile = $configHash->{'ORFFastaFileName'};#fasta file for the peptidic sequences created by dna2pep
 $outFolder = $configHash->{'outFolder'};#folder with the output files
 $statsFolder = $configHash->{'statsFolder'};#folder with statistics
 my $sessionFolder = $configHash->{'sessionFolder'};#folder with all the files for a given analysis
 $logFolder = $sessionFolder."/".$configHash->{'logFolder'};#folder with all the logs
 
 my $NCThresh = $configHash->{'NCThresh'};#Threshold of the potential score to be non coding
 my $NCORFLength = $configHash->{'NCORFLength'};#Maximum lenght of the longest ORF to be a non coding sequence
 my $NCSeqLength = $configHash->{'NCSeqLength'};#The sequence must be at least this lenght to be non coding

 my $dbDataFolder = $configHash->{'dbDataFolder'};#Folder with the database files
 my $mySqlUser = $configHash->{'mySqlUser'};#MySQLdb user id
 my $mySqlPass = $configHash->{'mySqlPass'};#MySQL db password
 $database = $configHash->{'database'};#MySQL database name
 my $outStatsFile;#Statistics file.
 
 $cdName4Expression = $configHash->{'cdName4Expression'};#The cd database to be used to plot the cd expression
 
 $separator = $configHash->{'separator'};#The separator used
 my $cdTableHash;#An hash to store informations about CD from file
 
 #CREATING connection to Uniprot DB
 my $platform = $configHash->{'platform'};
 my $host = $configHash->{'host'};## MySQL Host
 my $port = $configHash->{'port'};
 my $dbargs = {AutoCommit => 0, PrintError => 1};
 $dsn = "dbi:$platform:$database:$host:$port";
 
 
  #Variables have to be set.. 
 $outFolder =  $sessionFolder.'/'. $outFolder;#building output folder path
 $statsFolder =  $sessionFolder.'/'. $statsFolder; #building stats folder path 
 $querySeqs = $sessionFolder.'/'.$querySeqs;#File of the sequences
 $outStatsFile = "stats_".$annocriptOutFile;#Statistics file. Don't remove from this position!!!
 $annocriptOutFile = $outFolder.'/'.$annocriptOutFile; #Annocript output
 
 $maxLengthDescs = $configHash->{'maxLengthDescs'};
 
 my $NCCount = 0;#Number of non-coding sequences
 
 my $pathwaysIdsHash;#An hash containing IDs and descriptions
 
 # Number of results taken for the expression. We used this to check if something is going wrong with
 # the decision to take a different database for Conserved Domains when we are collecting abundances
 # for expression
  my $expResTaken = 0;

  #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;

 #Opening the Annocript output file. With this file we will extract statistics
 print "\n Now reading the Annocript output table for write statistics...\n";
 open(IN_FILE,"<$annocriptOutFile") or die "ERROR [$!]: Unable to open $annocriptOutFile. Please restart Annocript with doBuildOutput = YES\n "; 

 
 #Configuring the variables to access the final table
 setFieldNames($outFolder."/".$configHash->{'outFileName'},$configHash->{'dbInUse'});
 #print "Try to print something from headHash: ".$headHash->{'TranscriptName'}."\n";
 #print "and now with the variable: ".$headHash->{$seq}."\n";
 
	#Only if the output of lncRNAprediction program is present then we can create a file with the non coding sequences
	#Look at lncRNA_fasta_file_generator subroutine to see what are the constrains to be non coding
	$configHash->{'lncRNAPredOutPresent'} = checkOutputsPresence($configHash->{'lncRNAPredOut'},$outFolder);
	if ( $configHash->{'lncRNAPredOutPresent'} eq 'YES'){
		print "Creating a file with the long non-coding RNA sequences and one for the coding selected by Annocript....";
		#$configHash->{'NCCount'} = lncRNA_fasta_file_generator($configHash->{'outFileName'},$configHash->{'NCOutFile'},$outFolder);
		$configHash->{'NCCount'} = lncRNA_coding_fasta_file_generator($configHash->{'outFileName'},$configHash->{'NCOutFile'},$configHash->{'codingOutFile'},$outFolder,$configHash->{'NCORFLength'});
	
		$partDuration = time - $partTime;
		print "Finished writing NC information file in\t".$partDuration." seconds\n";
		$partTime = time;
	}

 #Create an hash with all the domains from the database
 createCDHash($dbDataFolder,$database,$mySqlUser,$mySqlPass,\$cdTableHash);
      
 #DB Connect
 my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
 if ($dbh->err()) { die "$DBI::errstr\n"; }
  
 
 #Variables for an array that will contain all the data
 my $totRows = 0;
 my $totCols = 0;
 my @outAnnArray = (); #Array with all Annocript values

 #Count annotated sequences
 my $totAnnotated  = 0;
 my $totFields = $headHash->{$NC4Annocript};
 #my $totFields = $seq;#If we are taking from the complete output
 my $totTranscripts = 0;
 
 
 #Here we store the output of Annocript in an array
 my $numField=0;#Counter of fields
 while (my $row = <IN_FILE>) {
  if ($totRows >0){
   my @fields = split(/\t/, $row);
   #print_array(\@fields);
   #print "Number of fields in row $totRows = ".scalar(@fields)."\n";#DEBUGCODE
   $numField=0;
   #Here we fill the array and count the fields
   foreach my $field (@fields){
     $outAnnArray[$totRows-1][$numField] = $field;
     $numField++;
     #print "$field\t";
   }
   #print "\n";
   #print "Fields on row ".($totRows + 1)." are: $numField \n ";#DEBUGCODE
   if ( countAnnotated(\@fields,$configHash->{'dbInUse'}) == 1 ){
        $totAnnotated++;
    }
    # numField is the right number while totFields starts from 0 
   if( ($numField-1) != $totFields){
      print "ERROR in output_and_stats: the number of fields found on the row ".($totRows + 1)." are $numField and not ".($totFields+1)." as expected.\n";
    }
  }
  $totRows++;
  
 }
 #print "final totRows = $totRows\n";#DEBUGCODE
 $totCols = $numField++;
 $totTranscripts = $totRows-1;#Number of transcripts are the number of rows without the header
 ########################################
 
 
 
  
 if ($totAnnotated == 0){
    print "WARNING in output_and_stats: something strange is happen. Please check why totAnnotated is zero...\n";
 }
   
 #A VERY LARGE LOOP TO EXTRACT INFORMATIONS
 #Variables needed
 my $maxSeqLength = -1;
 my $maxSeqLengthInd = -1;#Index of the max length sequence

 my $minSeqLength = 100000000;
 my $minSeqLengthInd = -1;#Index of the min length sequence

 my $meanSeqLength = 0;
 my $meanSeqLengthInd = -1;#Index of the mean length sequence

 my $medianSeqLength = 0;#You can order the lengths and take the totTranscripts/2 element
 my $medianSeqLengthInd = -1;#Index of the median length sequence


 
 #An hash will be used to sort the elements based on their occurrences
 my $bpHash;#biological process
 my $mfHash;#molecular function
 my $ccHash;#cellular component
 my $bpHashDom;#biological process for domains
 my $mfHashDom;#molecular function for domains
 my $ccHashDom;#cellular component for domains 
 my $CDHash;#domains
 my $PWHashL1;#level1 Pathways
 my $PWHashL2;#level2 pathways
 my $PWHashL3;#level3 pathways
 my $nearOrg1;#An hash for the nearest organism occurrences for SwissProt
 my $nearOrg2;#An hash for the nearest organism occurrences for UniRef/TrEMBL
	
 #the statistics will come out from the output file... and what happens if something is not there?
 #here we use some variables to be sure that each result is present.
 #If it is not there, the corresponding file will not be created
 my $trsWithBP = 0;
 my $trsWithMF = 0;
 my $trsWithCC = 0;
 #The same but for domains
 my $trsWithBPDom = 0;
 my $trsWithMFDom = 0;
 my $trsWithCCDom = 0;
 my $trsWithDom = 0;
 my $trsWithPwL1 = 0;
 my $trsWithPwL2 = 0;
 my $trsWithPwL3 = 0;
 
 my $existOrganisms1 = 0;
 my $existOrganisms2 = 0;
 
 #Check if the counts are present and save it in a variable
 #COUNTS ARE ACTUALLY NOT USED AND COMMENTED FOR FUTURE VERSION ADDITIONS
 my $countsPresent = 0;#If the counts are written in the output it will be '1'

 #$partDuration = time-$partTime;
 #print "Just finished with all the settings to create statistics in ".$partDuration." seconds\n\n";
 #$partTime = time;
      
 #Variables to take in count the strand information
 my $strand1Pos = 0;
 my $strand2Pos = 0;
 my $strand1Neg = 0;
 my $strand2Neg = 0;
 
 #Number of results counting
 my $totProteins1 = 0;
 my $totProteins2 = 0;
 my $totDomains = 0;
 my $totRibosomials = 0;
   
 print "Starting the loop to extract informations...\n";#DEBUGCODE
 #the loop can start
 my $iRow;
 for ($iRow = 0; $iRow < $totTranscripts; $iRow++ ){
  #Extract statistics
  
  #Maximum length
  if ( $outAnnArray[$iRow][$headHash->{$transcrLength}] > $maxSeqLength){
   $maxSeqLengthInd = $iRow;
   $maxSeqLength = $outAnnArray[$iRow][$headHash->{$transcrLength}];
  }
  #Minimum length
  if ( $outAnnArray[$iRow][$headHash->{$transcrLength}] < $minSeqLength){
   $minSeqLengthInd = $iRow;
   $minSeqLength = $outAnnArray[$iRow][$headHash->{$transcrLength}];
  }
  #Mean lenght of the sequences computation
  $meanSeqLength = $meanSeqLength + $outAnnArray[$iRow][$headHash->{$transcrLength}];
  
  #Compute GO stats for proteins  
  if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for biological processes
	  if($outAnnArray[$iRow][$headHash->{$BPId}] ne "-"){  
	    $trsWithBP++;#If at least one result is present this will be increased
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$BPId},\$bpHash,$countsPresent,undef);
	  }
	  
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for molecular functions
	  if($outAnnArray[$iRow][$headHash->{$MFId}] ne "-"){
	    $trsWithMF++;
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$MFId},\$mfHash,$countsPresent,undef);
	  }
	  
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for cellular components
	  if($outAnnArray[$iRow][$headHash->{$CCId}] ne "-"){
	    $trsWithCC++;
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$CCId},\$ccHash,$countsPresent,undef);
	  }
  }
  
  #Compute GO stats for domains
  if ( ($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both')){
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for biological processes
	  if($outAnnArray[$iRow][$headHash->{$domBPId}] ne "-"){  
	    $trsWithBPDom++;#If at least one result is present this will be increased
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$domBPId},\$bpHashDom,$countsPresent,'domains');
	  }
	  
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for molecular functions
	  if($outAnnArray[$iRow][$headHash->{$domMFId}] ne "-"){
	    $trsWithMFDom++;
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$domMFId},\$mfHashDom,$countsPresent,'domains');
	  }
	  
	  #Extract GO terms informations and put them in an hash
	  #It uses the array built before. Here for cellular components
	  if($outAnnArray[$iRow][$headHash->{$domCCId}] ne "-"){
	    $trsWithCCDom++;
	    extractGOstats($dbh,\@outAnnArray,$iRow,$headHash->{$domCCId},\$ccHashDom,$countsPresent,'domains');
	  }  
  }
  #print "gofinished: $iRow\n";
  # Here domains informations are written in the hash table
  # We decided to keep as description only ten words of the total 
  # description (very long) of the domain. The setting is changeable here with NumWords
  if($outAnnArray[$iRow][$headHash->{$CDName}] ne "-"){
   $trsWithDom++;
   extractCDstats($dbh,\@outAnnArray,$iRow,\$CDHash,\$expResTaken,$countsPresent);
  }

  # Here pathways informations are written in the hash table
  # We decided to keep as description only ten words of the total 
  # description (very long) of the pathways. The setting is changeable with the variable maxLengthDescs in the config file
  #if($outAnnArray[$iRow][$headHash->{$PwLev1}] ne "-"){
  if( ($outAnnArray[$iRow][$headHash->{$PwLev1}] !~ /^[\[\-\]]+$/) and ($outAnnArray[$iRow][$headHash->{$PwLev1}] ne "-") ){
   $trsWithPwL1++;
   extractPWstatsFromCols(\@outAnnArray,$iRow,$headHash->{$PwLev1},\$PWHashL1,$countsPresent);
  }
  # Here pathways informations are written in the hash table
  # We decided to keep as description only ten words of the total 
  # description (very long) of the pathways. The setting is changeable here with the variable maxLengthDescs in the config file
  #if($outAnnArray[$iRow][$headHash->{$PwLev2}] ne "-"){
  if( ($outAnnArray[$iRow][$headHash->{$PwLev2}] !~ /^[\[\-\]]+$/) and ($outAnnArray[$iRow][$headHash->{$PwLev2}] ne "-") ){
   $trsWithPwL2++;
   extractPWstatsFromCols(\@outAnnArray,$iRow,$headHash->{$PwLev2},\$PWHashL2,$countsPresent);
  }
  # Here pathways informations are written in the hash table
  # We decided to keep as description only ten words of the total 
  # description (very long) of the pathways. The setting is changeable here with the variable maxLengthDescs in the config file
  #if($outAnnArray[$iRow][$headHash->{$PwLev3}] ne "-"){
	if( ($outAnnArray[$iRow][$headHash->{$PwLev3}] !~ /^[\[\-\]]+$/) and ($outAnnArray[$iRow][$headHash->{$PwLev3}] ne "-") ){
   $trsWithPwL3++;
   extractPWstatsFromCols(\@outAnnArray,$iRow,$headHash->{$PwLev3},\$PWHashL3,$countsPresent);
  }

  #print "extractCDstats gofinished: $iRow\n";
  #Extracting an hash with the abundances of organisms from SwissProt
  if($outAnnArray[$iRow][$headHash->{$closerOS1}] ne "-"){
     $existOrganisms1++;
     my $organism = $outAnnArray[$iRow][$headHash->{$closerOS1}] ;
     $nearOrg1->{$organism}++;
  }

  #Extracting an hash with the abundances of taxonomies from UniRef or organisms from TrEMBL
  if($outAnnArray[$iRow][$headHash->{$closerOS2}] ne "-"){
     $existOrganisms2++;
     my $organism = $outAnnArray[$iRow][$headHash->{$closerOS2}] ;
     $nearOrg2->{$organism}++;
  }
    
  #Counting the number of non coding sequences from the output given
  if (defined $outAnnArray[$iRow][$headHash->{$NC4Annocript}]){
     if ( $outAnnArray[$iRow][$headHash->{$NC4Annocript}] == 1 ){
       $NCCount++;
     }
  }
  
  ########Here I count number of results
  #Counting the number of results from SwissProt
  if ( $outAnnArray[$iRow][$headHash->{$Name1}] ne "-"){
		 $totProteins1++;
  }
  
  #Counting the number of results from UniRef/TrEMBL
  if ( $outAnnArray[$iRow][$headHash->{$Name2}] ne "-"){
		 $totProteins2++;
	 }
  
  #Counting the number of results from Domains
  if ( $outAnnArray[$iRow][$headHash->{$CDName}] ne "-"){
		 $totDomains++;
	 }
  
  #Counting the number of results from Ribosomial
  if ( $outAnnArray[$iRow][$headHash->{$riboName}] ne "-"){
		 $totRibosomials++;
	 }
  
  #Counting the number of positive strands  from the first output given
  if ( $outAnnArray[$iRow][$headHash->{$Strand1}] ne "-"){
     if ( $outAnnArray[$iRow][$headHash->{$Strand1}] == 1 ){
			$strand1Pos++;
     }elsif( $outAnnArray[$iRow][$headHash->{$Strand1}] == -1 ){
			$strand1Neg++;
		 }
  }
  
  #Counting the number of positive strands  from the first output given
  if ( $outAnnArray[$iRow][$headHash->{$Strand2}] ne "-"){
     if ( $outAnnArray[$iRow][$headHash->{$Strand2}] eq "1" ){
				$strand2Pos++;
     }elsif( $outAnnArray[$iRow][$headHash->{$Strand2}] eq "-1" ){
				$strand2Neg++;
		}
  }
  
 }
  #print "extract organisms gofinished: $iRow\n";
 #Close db connection
  $dbh->commit();
  $dbh->disconnect();
 #########################CHECKING ERRORS###########################
  if( ($countsPresent == 1) and ($expResTaken == 0) ){
    print "WARNING!! Pay attention: in the plot for Conserved domains expression no result will compare. This means that"
        ." perhaps you chose, as cdName4Expression in the Annocript configuration, something that is not existent in the Cdd "
        ." (conserved domains database).\n";
  }
 $partDuration = time-$partTime;
 print "Finished with the creation of the statistics in\t".$partDuration." seconds\n\n";
 $partTime = time;
  
 #########################STATISTICS PRINTING########################
 
 #Here we print various statistics on a text file. Then this file will be used also to create an HTML page.
 print "Printing stats on a file.. it will need a while...\n";
 #Opening a file that will contain all the statistics
 open(OUT_STATS,">$outStatsFile");

 #The mean sequences length
 $meanSeqLength = $meanSeqLength/$totTranscripts;
 
 print OUT_STATS "The file of sequences is $querySeqs\n";
 print OUT_STATS "The total number of sequences is $totTranscripts\n";
 print OUT_STATS "The mean sequence length is ".int($meanSeqLength)."\n";

 #Minimum and maximum sequences length
 print OUT_STATS "The minimum and maximum sequence length are respectively $outAnnArray[$minSeqLengthInd][$headHash->{$transcrLength}]". 
    " and $outAnnArray[$maxSeqLengthInd][$headHash->{$transcrLength}]\n";
 
 #Here we set an header to write in the GO file table
 my @headerPerc = ('rowind','go_id','occurrences','tot_count','description');#The header
 my @headerExpr = ('rowind','go_id','occurrences','description');#The header
 my $dataName = '';#Name of the data that we are going to print
 my $tableName = '';#Name of the table that we are going to print
 
 #
 if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){ 
	 #Print the go bp counts on a file 
	 $dataName = 'biological processes';
	 if ($trsWithBP>0){
	   $tableName = $configHash->{'bpStatsFile'};
	   printOutGOStats(\@headerPerc,\@headerExpr,$bpHash,$tableName,$trsWithBP,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	   #print Dumper($bpHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
	 
	 #Print the go mf counts on a file
	  $dataName = 'molecular functions'; 
	 if ($trsWithMF>0){
	   $tableName = $configHash->{'mfStatsFile'};
	    printOutGOStats(\@headerPerc,\@headerExpr,$mfHash,$tableName,$trsWithMF,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	    #print Dumper($mfHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
	 
	 #Print the go cc counts on a file
	 $dataName = 'cellular components';
	 if ($trsWithCC>0){
	   $tableName = $configHash->{'ccStatsFile'};
	   printOutGOStats(\@headerPerc,\@headerExpr,$ccHash,$tableName,$trsWithCC,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	   #print Dumper($ccHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
 }
 
 if ( ($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both')){
	 #Print the go bp counts on a file 
	 $dataName = 'biological processes';
	 if ($trsWithBPDom>0){
	   $tableName = $configHash->{'bpStatsFileDom'};
	   printOutGOStats(\@headerPerc,\@headerExpr,$bpHashDom,$tableName,$trsWithBPDom,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	   #print Dumper($bpHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
	 
	 #Print the go mf counts on a file
	  $dataName = 'molecular functions'; 
	 if ($trsWithMFDom>0){
	   $tableName = $configHash->{'mfStatsFileDom'};
	    printOutGOStats(\@headerPerc,\@headerExpr,$mfHashDom,$tableName,$trsWithMFDom,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	    #print Dumper($mfHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
	 
	 #Print the go cc counts on a file
	 $dataName = 'cellular components';
	 if ($trsWithCCDom>0){
	   $tableName = $configHash->{'ccStatsFileDom'};
	   printOutGOStats(\@headerPerc,\@headerExpr,$ccHashDom,$tableName,$trsWithCCDom,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
	   #print Dumper($ccHash);
	 }else{print "There are no $dataName classes in the output file. The statistics will not be created.\n";}
 }
#Here we set an header to write in the pathways file table
 @headerPerc = ('rowind','pw_id','occurrences','tot_count','description');#The header
 @headerExpr = ('rowind','pw_id','occurrences','description');#The header
 
 #Print the pathways first level on a file
  $dataName = 'pathways first level';
 if ($trsWithPwL1 > 0){
   $tableName = $configHash->{'pathwaysL1File'};
   printOutGOStats(\@headerPerc,\@headerExpr,$PWHashL1,$tableName,$trsWithPwL1,$statsFolder,\@outAnnArray,$dataName,$countsPresent);
   #print Dumper($PWHashL1);
 }else{print "There are no $dataName in the output file. The statistics will not be created.\n";}


 #Print the pathways first level on a file
  $dataName = 'pathways second level';
 if ($trsWithPwL2 > 0){
   $tableName = $configHash->{'pathwaysL2File'};
   printOutGOStats(\@headerPerc,\@headerExpr,$PWHashL2,$tableName,$trsWithPwL2,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
   #print Dumper($PWHashL2);
 }else{print "There are no $dataName in the output file. The statistics will not be created.\n";}


 #Print the pathways first level on a file
  $dataName = 'pathways third level';
 if ($trsWithPwL3 > 0){
   $tableName = $configHash->{'pathwaysL3File'};
   printOutGOStats(\@headerPerc,\@headerExpr,$PWHashL3,$tableName,$trsWithPwL3,$statsFolder,\@outAnnArray,$dataName,$countsPresent); 
   #print Dumper($PWHashL3);
 }else{print "There are no $dataName in the output file. The statistics will not be created.\n";}

  #Print the domains counts on a file
  if ($trsWithDom>0){
    #Here we set an header to write in the GO file table
    my @headerPerc = ('rowind','go_id','occurrences','tot_count','description');#The header
    my @headerExpr = ('rowind','go_id','occurrences','description');#The header
    $dataName = 'domains';
     #Create percentages respect to the number of Annotated sequences
    foreach my $key (%$CDHash){
      if (defined $CDHash->{$key}){
       $CDHash->{$key}->{'perc'} = (($CDHash->{$key}->{'count'})/$trsWithDom)*100;
      }
    }
    #print Dumper($CDHash);
    @headerPerc = ('rowind','cd_id','occurrences','tot_count','description');#The header
     my $filePerc = $configHash->{'cdStatsFile'}."_percentages";
    sortAndWrite($filePerc,$CDHash,\@headerPerc,"perc");
    print "A file for $dataName percentages has been created...\n";
    move($filePerc,$statsFolder) or print "ERROR [$!]: unable to move $filePerc in $statsFolder. Please check permissions.\n";

  }else{print "There are no domains info in the output file. The statistics will not be created.\n";}
  
 #Print the organisms counts on a file
 if ($existOrganisms1 > 0){
  my @header;
  @header = ('row_ind','description','occurrences');#The header
  sortAndWriteSingle($configHash->{'closerOSFile1'},$nearOrg1,\@header);
  move($configHash->{'closerOSFile1'},$statsFolder) or print "ERROR [$!]: unable to move ".$configHash->{'closerOSFile1'}." in $statsFolder. Please check permissions.\n";
 }else{print "There are no organisms info in the output file. The statistics will not be created.\n";}
 
 #Print the organisms counts on a file
 if ($existOrganisms2 > 0){
  my @header;
  @header = ('row_ind','description','occurrences');#The header
  sortAndWriteSingle($configHash->{'closerOSFile2'},$nearOrg2,\@header);
  move($configHash->{'closerOSFile2'},$statsFolder) or print "ERROR [$!]: unable to move ".$configHash->{'closerOSFile2'}." in $statsFolder. Please check permissions.\n";
 }else{print "There are no organisms info in the output file. The statistics will not be created.\n";}
 
 
 #This stats are done only if the alphabet is with ATGC chars
 if ($configHash->{'sequencesAlphabet'} eq 'dna'){
   #We extract percentages of presences of ATGC in the fasta file
   my @percentages = extractPercentages($querySeqs,$totTranscripts,"ATGC_Percentages.txt");
   print "A file for percentages of ATGC has been created...\n";
   move($configHash->{'ATGCPercFile'},$statsFolder) or print "ERROR [$!]: unable to move ".$configHash->{'ATGCPercFile'}." in $statsFolder. Please check permissions.\n";
   
   printf OUT_STATS "Average percentage of Adenine: %.2f \n Average percentage of Guanine: %.2f \n",$percentages[0],$percentages[1];
   printf OUT_STATS "Average percentage of Thymine: %.2f \n Average percentage of Cytosine: %.2f \n",$percentages[2],$percentages[3];
   printf OUT_STATS "Average percentage of N: %.2f \n Average percentage of GC: %.2f \n ",$percentages[4],$percentages[5];
 }else{
   print "The sequences are not in DNA format so Annocript will not perform the statistics about ATGC percentages\n";
   }
 
 
  if( check_presence( $outFolder."/". $configHash->{'dna2pepOut'}) ){
     $configHash->{'dna2pepOutPresent'} = 'YES';
  }else{
    $configHash->{'dna2pepOutPresent'} = 'NO';
  }
  
  my $secondDB = 'UniRef';
  if ($configHash->{'dbInUse'} ne 'uniref'){
		$secondDB = 'TrEMBL';
  }
  
  #Number of results from the blast
  print OUT_STATS "Number of blast results:\n";
  if ($totProteins1 > 0) {print OUT_STATS "Swiss-Prot: $totProteins1\n";}
  if ($totProteins2 > 0) {print OUT_STATS "$secondDB : $totProteins2\n";}
  if ($totDomains > 0) {print OUT_STATS "Domains: $totDomains\n";}
  if ($totProteins1 > 0) {print OUT_STATS "Ribosomal RNAs: $totRibosomials\n";}
 print OUT_STATS "\n";
 
  #Information about strand, if present
  if (($strand1Pos > 0) or ($strand2Pos > 0)) {
		 print OUT_STATS "Alignments in the positive strand:\n";
	} 
  if ($strand1Pos > 0) {print OUT_STATS "Swiss-Prot: $strand1Pos\n";}
  if ($strand2Pos > 0) {print OUT_STATS "$secondDB: $strand2Pos\n";}
  print OUT_STATS "\n";
  
  if (($strand1Neg > 0) or ($strand2Neg > 0)) {
		print OUT_STATS "Alignments in the negative strand:\n";
	}
  if ($strand1Neg > 0) {print OUT_STATS "Swiss-Prot: $strand1Neg\n";}
  if ($strand2Neg > 0) {print OUT_STATS "$secondDB: $strand2Neg\n";}
  print OUT_STATS "\n";
	print OUT_STATS "Transcripts with at least one blast result: ".$totAnnotated."\n\n";
   
  #Detection of the Agreement: to see if the strand of the longest ORF is the same of the one in blast
  if ( $configHash->{'dna2pepOutPresent'} eq 'YES'){
    print OUT_STATS "Results in agreement with strand of the longest ORF: ".detectLORFAgreement($orfFile,$outFolder)."\n";
  }
 print OUT_STATS "\n";
 
 #Printing the number of non coding sequences
 if( check_presence( $outFolder."/".$configHash->{'lncRNAPredOut'})){
     $configHash->{'lncRNAPredOutPresent'} = 'YES';
  }
  if ($NCCount > 0){
    print OUT_STATS "Number of non coding sequences: $NCCount \n (obtained with probability lower than: $NCThresh and ".
    "maximum length of the orf: $NCORFLength)\n\n";
  }

 #Close the file of statistics. It will be used later to create an html page
 close (OUT_STATS);
 
 #Move the final stats file to the stats folder
 move($outStatsFile,$statsFolder);
 print "File $outStatsFile correctly written in $statsFolder...\n";
 
 # An html file with plots and data is printed.
 plots_2_html($sessionFolder,$outStatsFile);
  
  $partDuration = time-$partTime;
  print "Finished writing an html file with statistics in\t".$partDuration." seconds.";
  $partTime = time;
  print " It has been stored in $outStatsFile\n\n";
 
  #Calculate the time needed for all the processes
	my $globalDuration = time - $globalStart;
	print "\n".scalar(localtime)." Annocript finished the creation of statistics in: ".$globalDuration." seconds\n";
}



=head2 extractCDstats

 Title : extractCDstats
 Usage : extractCDstats( -dbh => 'Reference to the database accessed',
                         -outAnnArray => 'the output of Annocript in a matrix',
                         -iRow => 'the current row analyzed',
                         -CDName => 'Number of the field in the table of the go ID'
                         -CDHash => 'The CD hash that we are building with the statistics'
                         
     );

 Function: this function builds an hash counting the occurrences of domains names

 Returns : nothing but the $CDHash inserted will be modified

=cut
sub extractCDstats {
  #Input
  my $dbh = shift;#Reference to the database accessed
  my $outAnnArray = shift;#The matrix with all the table of output of Annocript
  my $iRow = shift;#the current row analyzed
  my ($CDHash) = shift;#The CD hash that we are building with the statistics
  my $expResTaken = shift;
  my $countsPresent = shift;
  
  #Local variables
   my $splitInd = 0;
   my @cdNamesFound = ();#Elements in the row that have been already inserted
   my @splitted = split(/\Q$separator\E/,${$outAnnArray}[$iRow][$headHash->{$CDName}]);#Array to split domain names
   
   #The loop goes for each of the elements separated
   foreach my $cdName (@splitted){
     if( $cdName =~ /$cdName4Expression/){
       #and if the element hasn't been already evaluated
      if( ! (any{/$cdName/} @cdNamesFound) ){
        #push it in the already examined list and access db to fetch its informations
       push(@cdNamesFound ,$cdName);
       my @res = ();
       my $select_handle;
       #Extracts information from UNIPROT only if $database is defined. (blastx i.e.) 
       if (defined $database){
          #fetch informations from the database
          $select_handle = $dbh->prepare ("SELECT cdName, cdShortName, cdDescription FROM cdTable WHERE cdName='".$cdName."'");
          die "ERROR [$!]: can't prepare queries; aborting" unless defined $select_handle;
          $select_handle->execute() or die $select_handle->errstr;
          @res = $select_handle->fetchrow_array;
       }else{print "ERROR [$!]: $database connession fault in a point where Annocript don't expect\n";}
       
       #If results are coming out from the database...
       if ( scalar(@res) > 0 ){
          my $desc = $res[0]." ".$res[1]." ".$res[2];
          
          ##..we take the description from the result and separate it to get only $numWords words.
           ##this is done to not have a description too long (it will be printed on an image with R)
          my $reducedDesc = substr $desc, 0 , $maxLengthDescs;
          
          $$CDHash->{$cdName}->{'count'}++;
          $$CDHash->{$cdName}->{'desc'} = $reducedDesc;
          
       }else{print "ERROR [$!]: for the domain: $cdName has not been found a corresponding domain in the database.\n";}
      }
   }
  }
}

=head2 extractPWstatsFromCols

 Title : extractPWstatsFromCols
 Usage : extractPWstatsFromCols( -dbh => 'Reference to the database accessed',
                         -outAnnArray => 'the output of Annocript in a matrix',
                         -iRow => 'the current row analyzed',
                         -PWName => 'Number of the field in the table of the go ID'
                         -PWHash => 'The PW hash that we are building with the statistics'
                         
     );

 Function: this function builds an hash counting the occurrences of pathways names

 Returns : nothing but the $PWHash inserted will be modified

=cut
sub extractPWstatsFromCols {
  #Input for the function
  my $outAnnArray = shift;#The matrix with all the table of output of Annocript
  my $iRow = shift;#the current row analyzed
  my $PwId = shift;#The index in the output of Annocript where the level is found
  my ($PWHash) = shift;#The hash that we are building with the statistics
  my $countsPresent = shift;#says if the counts are present to get the enrichments
  
  #Local variables  
  my @pwFound = ();#Elements in the field that have been already inserted

  my @splitted = split(/\Q$separator\E/,${$outAnnArray}[$iRow][$PwId]);#Array to split different pathways associated   
     
   #The loop goes for each of the elements separated
   foreach my $pathway (@splitted){
     #if it is not there
     if( (!( any {/\Q$pathway\E/} @pwFound)) and ($pathway ne "-") ){
        #push it in the already examined list and access hash to fetch its informations
        push(@pwFound ,$pathway);
            
        my $reducedDesc = substr $pathway, 0 , $maxLengthDescs;
        $$PWHash->{$pathway}->{'count'}++;
        $$PWHash->{$pathway}->{'desc'} = $reducedDesc;
     }
  }
}



=head2 printOutGOStats

 Title : printOutGOStats
 Usage : printOutGOStats(  -header => 'the header to write in the file',
                           -GOHash => 'the hash where the info are stored',
                           -fileName => 'name of the file to print',
                            -totAnnotated => 'number of annotated sequences',
                            - statsFolder => 'path to the stats folder',
                            - outAnnArray => 'array storing the annocript out table'
                            - dataName => 'data name to print'
            );

 Function: Create percentages respect to the number of Annotated sequences and write them sorted in a given 
            output file

 Returns : nothing

=cut  
sub printOutGOStats {
  my $headerPerc = shift;
  my $headerExpr = shift;
  my $GOHash = shift;
  my $fileName = shift;
  my $totAnnotated = shift;
  my $statsFolder = shift;
  my $outAnnArray = shift;
  my $dataName = shift;
  my $countsPresent = shift;
  
  
    #This counts the simple abundances and is commented
    #sortAndWrite($fileName,$GOHash,$header,"count");
    #print "A file for biological process has been created...\n";
    #move($fileName,$statsFolder) or print "ERROR [$!]: unable to move $fileName in $statsFolder. Please check permissions.\n";
    
    #Create percentages respect to the number of Annotated sequences
    foreach my $key (%$GOHash){
      if (defined $GOHash->{$key}){
       $GOHash->{$key}->{'perc'} = (($GOHash->{$key}->{'count'})/$totAnnotated)*100;
      }
    }
    
    #print Dumper($bpHash);
    my $filePerc = $fileName."_percentages";
    sortAndWrite($filePerc,$GOHash,$headerPerc,"perc");
    print "A file $fileName with $dataName percentages has been created...\n";
    move($filePerc,$statsFolder) or print "ERROR [$!]: unable to move $filePerc in $statsFolder. Please check permissions.\n";
    
}


=head2 plots_2_html

 Title : plots_2_html
 Usage : plots_2_html(  -sessionFolder => 'the folder where to find the stats one',
                        -outStatsFile => 'the folder where the stats are kept'
            );

 Function: Goes through all the files in the folder stats and matches the names
					with specific words. Using these words choses the names and uses R to 
					print on JPG images a pie chart or a barplot.

 Returns : nothing

=cut  
sub plots_2_html{
 my $sessionFolder = shift;
 my $outStatsFile = shift;
 
 my $statsFolder = $sessionFolder.'/stats';#Statistics folder
 $outStatsFile = $statsFolder."/".$outStatsFile;#Statistics file
 
 #Setting paths for the folder and the html page
 my $imgDir = $configHash->{'imgDir'};#The name of the directory (to use for links)
 my $htmlHome =  $configHash->{'htmlHome'};
 my $imgPath = $configHash->{'imgPath'};#the path to the directory (to access for analysis)
 
 #Number of elements to show in the plots
 my $topToShow = $configHash->{'topToShow'};
 
 #Strings to match to identify files
 my $cc = $configHash->{'ccStatsFile'}.'_';
 my $mf = $configHash->{'mfStatsFile'}.'_';
 my $bp = $configHash->{'bpStatsFile'}.'_';
 my $ccDom = $configHash->{'ccStatsFileDom'}.'_';
 my $mfDom = $configHash->{'mfStatsFileDom'}.'_';
 my $bpDom = $configHash->{'bpStatsFileDom'}.'_';
 my $CD = $configHash->{'cdStatsFile'}.'_';
 my $os1 = $configHash->{'closerOSFile1'};
 my $os2 = $configHash->{'closerOSFile2'}; 
 my $pwl1 = $configHash->{'pathwaysL1File'}.'_';
 my $pwl2 = $configHash->{'pathwaysL2File'}.'_';
 my $pwl3 = $configHash->{'pathwaysL3File'}.'_';
 
 my $table = 'table';
 my $expression = 'expression';
 my $percentages = 'percentages';
 
 #Script in R to create the plots
 my $R_plot_script;
 #Here the script for the plots will be changed according to the user choice
 #if ($configHash->{'plotType'} eq 'piechart'){
	# $R_plot_script = $configHash->{'ProgOutParserFolder'}."/".$configHash->{'R_piechart_script'};
 #}else{
	 $R_plot_script = $configHash->{'ProgOutParserFolder'}."/".$configHash->{'R_barplot_script'};
	 my $plotType = "barplot";
 #}
  
 my $R_seqs_stats_script = $configHash->{'ProgOutParserFolder'}."/".$configHash->{'R_seqs_stats_script'};
 
 my $RLogPath = $logFolder."/".$configHash->{'R_log_file'};

 #Creating the images dir
 if(!(-e $imgPath )){
  mkdir $imgPath or die "ERROR [$!]: cannot create dir $imgPath\n";
 }else {print "$imgPath already exists...\n";}
 
 #Move the logo inside the images dir
 copy($configHash->{'usefulFolder'}."/".$configHash->{'AnnocriptLogo'},$imgPath."/".$configHash->{'AnnocriptLogo'}) 
	or print "Cannot copy the logo in the $imgPath folder\n";
 
 
 #
 opendir DIR, $statsFolder or die "ERROR [$!]: cannot open dir $statsFolder\n";
 my @files= readdir DIR;
 closedir DIR;
 my $workingDir = getcwd;

	my $dbChosen = '';
 if ($configHash->{'dbInUse'} eq 'uniref'){
		$dbChosen = 'UniRef';
	}else{
		$dbChosen = 'TrEMBL';	
	}
  #For each file in the images folder if its name contains the word 'table'
  #then it is a table that has to be plot
 foreach my $file (@files){
	 #print "Printing a plot in R with:  $file...\n";#DEBUGCODE
  if ($file =~ /$table/){
   #chdir $dir;
   
   my $name;
   #IMPORTANT: Names should be one single word or the R command will fail!!
   if ( $file =~ /$cc/ ){
    $name = 'Cellular_Component_abundances';
   }elsif ( $file =~ /$mf/ ){
    $name = 'Molecular_Function_abundances';
   }elsif ( $file =~ /$bp/ ){
    $name = 'Biological_Process_abundances';
   }elsif ( $file =~ /$ccDom/ ){
    $name = 'Cellular_Component_abundances_domains';
   }elsif ( $file =~ /$mfDom/ ){
    $name = 'Molecular_Function_abundances_domains';
   }elsif ( $file =~ /$bpDom/ ){
    $name = 'Biological_Process_abundances_domains';
   }elsif ( $file =~ /$pwl1/ ){
    $name = 'Level1_pathways_abundances';
   }elsif ( $file =~ /$pwl2/ ){
    $name = 'Level2_pathways_abundances';
   }elsif ( $file =~ /$pwl3/ ){
    $name = 'Level3_pathways_abundances';
   }elsif ( $file =~ /$CD/ ){
    $name = 'Conserved_Domains_abundances';
   }elsif ( $file =~ /$os1/ ){
    $name = 'Closer_organisms_abundances_(SwissProt)';
   }elsif ( $file =~ /$os2/ ){
    $name = "Closer_organisms_abundances_($dbChosen)";
   }else{
    $name = $plotType.":".$file;
   }
 
 
   my $outImg = $imgPath."/".$file.".jpg";
   my $xlab = 'GO_Term';
   my $ylab = 'Expression_level'; 
   $file = $statsFolder."/".$file;

   #BAR CHARTS ARE ACTIVE
   my $command = "R CMD BATCH --no-save --no-restore '--args $file $outImg $name $xlab $ylab $topToShow' $R_plot_script $RLogPath";
   #print "The command is:\n $command\n";
   ( system ($command) ) == 0
     #  or die "ERROR [$!]: an error occurred while running R with R CMD BATCH \n";
      or R_die($command,$R_plot_script);
   #chdir $workingDir;
  } 
  
  #If the file name contains the word enrichment than the name is different
  if ($file =~ /$expression/){
   #  chdir $dir;
   #print "Opening $file...\n";
   my $name;
   #IMPORTANT: Names should be one single word or the R command will fail!!
   if ( $file =~ /$cc/ ){
    $name = 'Cellular_Component_expression';
   }elsif ( $file =~ /$mf/ ){
    $name = 'Molecular_Function_expression';
   }elsif ( $file =~ /$bp/ ){
    $name = 'Biological_Process_expression';
   }elsif ( $file =~ /$ccDom/ ){
    $name = 'Cellular_Component_expression_domains';
   }elsif ( $file =~ /$mfDom/ ){
    $name = 'Molecular_Function_expression_domains';
   }elsif ( $file =~ /$bpDom/ ){
    $name = 'Biological_Process_expression_domains';
   }elsif ( $file =~ /$pwl1/ ){
    $name = 'Level1_pathways_expression';
   }elsif ( $file =~ /$pwl2/ ){
    $name = 'Level2_pathways_expression';
   }elsif ( $file =~ /$pwl3/ ){
    $name = 'Level3_pathways_expression';
   }elsif ( $file =~ /$CD/ ){
    $name = 'Conserved_Domains_expression';
   }elsif ( $file =~ /$os1/ ){
    $name = 'Closer_organisms_expression_(SwissProt)';
   }elsif ( $file =~ /$os2/ ){
    $name = "Closer_organisms_expression_($dbChosen)";
   }else{
    $name = $plotType.":".$file;
   }
  
   #Printing the images with R
   my $outImg = $imgPath."/".$file.".jpg";
   $file = $statsFolder."/".$file;
   my $xlab = 'GO_Term';
   my $ylab = 'Percentages_of_abundance';
   #my $outImg = $file.".jpg";
   my $command = "R CMD BATCH --no-save --no-restore '--args $file $outImg $name $xlab $ylab $topToShow' $R_plot_script $RLogPath";
   #print "The command is:\n $command\n";
   ( system ($command) ) == 0
      # or die "ERROR: an error occurred while running R with R CMD BATCH ";
       or R_die($command,$R_plot_script);
   #chdir $workingDir;
  } 
  
  #If name  contains the word 'percentages' the name is different
  if ($file =~ /$percentages/){
   #  chdir $dir;
   #print "Opening $file...\n";
   my $name;
   #IMPORTANT: Names should be one single word or the R command will fail!!
   if ( $file =~ /$cc/ ){
    $name = 'Cellular_Component_percentages';
   }elsif ( $file =~ /$mf/ ){
    $name = 'Molecular_Function_percentages';
   }elsif ( $file =~ /$bp/ ){
    $name = 'Biological_Process_percentages';
   }elsif ( $file =~ /$ccDom/ ){
    $name = 'Cellular_Component_percentages_domains';
   }elsif ( $file =~ /$mfDom/ ){
    $name = 'Molecular_Function_percentages_domains';
   }elsif ( $file =~ /$bpDom/ ){
    $name = 'Biological_Process_percentages_domains';
   }elsif ( $file =~ /$pwl1/ ){
    $name = 'Level1_pathways_percentages';
   }elsif ( $file =~ /$pwl2/ ){
    $name = 'Level2_pathways_percentages';
   }elsif ( $file =~ /$pwl3/ ){
    $name = 'Level3_pathways_percentages';
   }elsif ( $file =~ /$CD/ ){
    $name = 'Conserved_Domains_percentages';
   }elsif ( $file =~ /$os1/ ){
    $name = 'Closer_organisms_percentages_(SwissProt)';
   }elsif ( $file =~ /$os2/ ){
    $name = "Closer_organisms_percentages_($dbChosen)";
   }else{
    $name = $plotType.":".$file;
   }
   my $outImg = $imgPath."/".$file.".jpg";
   $file = $statsFolder."/".$file;
   my $xlab = 'GO_Term';
   my $ylab = 'percentages_Amount';
   #my $outImg = $file.".jpg";
   
   #And finally we print with R
   my $command = "R CMD BATCH --no-save --no-restore '--args $file $outImg $name $xlab $ylab $topToShow' $R_plot_script $RLogPath";
   #print "The command is:\n $command\n";
   ( system ($command) ) == 0
      or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file in the log folder\n";
       #or R_die($command,$R_plot_script);
  }
 }
 
 #Generating plots with statistics
 my $annocriptOutFile = $outFolder.'/'.$configHash->{'outFiltered'};#Output file of Annocript
 if ( (-e $annocriptOutFile) and !(-z $annocriptOutFile) ){
	 my $name = $configHash->{'name'};#organism name
	 my $RlogFile2 = $logFolder."/R_seqs_stats_script.log";
	 my $command = "R CMD BATCH --no-save --no-restore '--args all $annocriptOutFile $imgPath $name' $R_seqs_stats_script $RlogFile2";
	 #print "The command is:\n $command\n";
		( system ($command) ) == 0
				or die "ERROR [$!]: an error occurred while executing $command. Please check the R log file in the log folder $logFolder\n";
	}

 #Print an HTML page with all the statistics about the transcriptome
 #Create a string with all the stats file content
 open FILE, $outStatsFile or die "ERROR [$!]: can't open file: $outStatsFile. Check permissions.\n"; 
 my $string = "<center>";
 while (<FILE>){
  $string = $string." <br>".$_;
 }
 close FILE;
 $string.="</center>";
 
 my @descText = ("Statistics for transcriptome","Homology statistics","Lengths and coverage");
 my @links = ($configHash->{'htmlHome'},$configHash->{'htmlAnnotStats'},$configHash->{'htmlSeqStats'});
 html_page($imgPath,\@descText,0,\@links,undef,$string);

 #HTML page with plots from homology results
 my $annotStatsSuffixes = $configHash->{'bpStatsFile'}."|".$configHash->{'mfStatsFile'}."|".$configHash->{'ccStatsFile'}."|".
 $configHash->{'bpStatsFileDom'}."|".$configHash->{'mfStatsFileDom'}."|".$configHash->{'ccStatsFileDom'}."|".
 $configHash->{'cdStatsFile'}."|".$configHash->{'closerOSFile1'}."|".$configHash->{'closerOSFile2'}."|".$configHash->{'pathwaysL1File'}."|".$configHash->{'pathwaysL2File'}."|".
 $configHash->{'pathwaysL3File'};
 $string = "";
 html_page($imgPath,\@descText,1,\@links,$annotStatsSuffixes,$string);
 
 #HTML page with coverage and lengths statistics
 my $imgsKeywords = "coverages.jpg|coverages.jpeg|lengths.jpg|lengths.jpeg";
 html_page($imgPath,\@descText,2,\@links,$imgsKeywords,$string);
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
			from ($descrText). Puts the logo of Annocript from the path. A loop searches for all the images in $imgDir
			with title comprising a word inside ($imgsKeywords). Print a series of links to the other pages using a loop
			searching in ($links) whose names are in descrText. Finally prints information at the bottom of the page.

 Returns : nothing

=cut 
sub html_page{
 my $imgPath = shift;#Directory where the images are stored	
 my ($descrText) = shift;#descriptive texts of the links to the other pages
 my $descrTextIndex = shift;#Index of the descriptive text of this page
 my ($links) = shift;#array with the links at the other pages
 my $imgsKeywords = shift;#keywords to use to search the images inside the folder
 my $infoText = shift;#A string to print as informations in the body of the page
 
  my $cgi = CGI->new; #Instantiate a CGI class
 
 #Page with Plots of length and coverages coming from the annotation
 open (OUT,">$statsFolder/".$$links[$descrTextIndex]) || die "ERROR [$!]: can't open file ".$$links[$descrTextIndex]."\n";
  
 print OUT $cgi->start_html("Annocript ".$configHash->{'AnnocriptVer'}."- ".$$descrText[$descrTextIndex]), # Begin HTML page
       #$cgi->header, # First we create a header
       $cgi->center; 
 print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$configHash->{'AnnocriptLogo'}."' ALT='image' ALIGN=TOP> <br>";
       #$cgi->h1('<i>Annocript '.$configHash->{'AnnocriptVer'}.'</i>'), # create the Tag <h1>TEXT</h1>
 print OUT $cgi->h3($$descrText[$descrTextIndex]."\n");
 print OUT $cgi->p($infoText);
 
 #Open images folder
 opendir IMG, $imgPath or die "ERROR [$!]: cannot open dir $imgPath.\n";
 my @images= readdir IMG;
 my @imgsSorted = sort {$a cmp $b} @images;#order alfabetically the images names
 closedir IMG;
 
 #Search and write HTML code to show images
 #print OUT "</center>";
 if(defined($imgsKeywords)){
	 foreach my $img (@imgsSorted){
	  #if (scalar(grep {/\b$img\b/} @$imgsKeywords) == 1){
	  if ($img =~ /$imgsKeywords/){
			print OUT "<IMG SRC='".extract_name($imgPath,0)."/".$img."' ALT='image' ALIGN=TOP> <br>";
	  } 
	 }
 }

 #print the website map
 for ( my $link=0; $link<scalar(@$links);$link++ ){
	print OUT $cgi->a({href=>$$links[$link]},$$descrText[$link]);
	print OUT " | ";
 }

 #Print information
 print OUT $cgi->h5('Annocript '.$configHash->{'AnnocriptVer'}.' - Copyright of Bioinformatics Lab SZN Naples'),
      $cgi->p(scalar localtime); # create <p>Text</p>
 print OUT "</center>";
 print OUT $cgi->end_html; # end of HTML
 close(OUT);
}


=head2 R_die

 Title : R_die
 Usage : R_die( );

 Function: die with R script and says to see the R log file. It moves the R log file in the log folder.

 Returns : nothing.. it dies..

=cut 
sub R_die{
  my $command = shift;
  my $R_plot_script = shift;
   
  print "Annocript will die for problems with R. Command: $command\n";
 
  if ( -d $logFolder){
   # print "Trying to move ".$R_plot_script."out in ".$logFolder."/".$R_plot_script."out\n";
    my $RlogName = extract_name($R_plot_script,0);
    #Moving the log file to the log folder before to exit Annocript
    #The working folder is Annocript so R will print the log file there...
    move($RlogName."out", $logFolder."/".$RlogName."out") 
      or die "Cannot move ".$R_plot_script."out in ".$logFolder."/".$R_plot_script."out";
    print "Please check the R log file at $logFolder/".$R_plot_script."out\n";
  }else{
    print "Annocript tried without succeed to move the log file ".$R_plot_script."out: $logFolder does not exists.\n";
  }
  die;
}  
  
=head2 countAnnotated

 Title : countAnnotated
 Usage : countAnnotated( -annocriptOutFile => 'the output of Annocript' );

 Function: this function counts the sequences with at least one result of blast programs from the Annocript's output

 Returns : the count of annotated sequences

=cut 
sub countAnnotated {
  my ($fields) = shift;
  my $databaseDl = shift;
  
  #print "fields[Name1]: ".$$fields[$headHash->{$Name1}]." fields[CDName] : ".$$fields[$headHash->{$CDName}]."\n";
 my $annotated = 0;
 
 if ($databaseDl eq 'uniprotkb'){
    #if there is at least one result from BLAST outputs...
    if ( !(($$fields[$headHash->{$Name1}] eq "-") and ($$fields[$headHash->{$Name2}] eq "-") and ($$fields[$headHash->{$CDName}] eq "-") and ($$fields[$headHash->{$riboName}] eq "-")) ){
     $annotated++;
    }
 }
 
 if ($databaseDl eq 'uniref'){
  #if there is at least one result from BLAST outputs...
    if ( !(($$fields[$headHash->{$Name1}] eq "-") and ($$fields[$headHash->{$Name2}] eq "-") and ($$fields[$headHash->{$CDName}] eq "-") and ($$fields[$headHash->{$riboName}] eq "-")) ){
     $annotated++;
    }
 }
 
 #print "Sequence annotated: $annotated\n";
 close (ANNOUT_FILE);
 return $annotated ;
}


=head2 sortAndWrite

 Title : sortAndWrite
 Usage : sortAndWrite( -fileName => 'the output ',
        -hash => 'an hash with the pairs GOID, occurrences',
         -header => 'the header to be written',
         field => ''
     );

 Function: this function starts from an hash with the informations about bp, mf, cc or domains and
   order it from the highest number of occurrences to the lowest. 
   All the pairs sorted are written in an input file

 Returns : a file with the sorted elements

=cut
sub sortAndWrite{
 my $fileName=shift;
 my $hash = shift;
 my $header =shift;
 my $field = shift;
 
 my $rowCount = 0;
 
 open(DESC, ">".$fileName);
 print DESC join("\t",@$header);
 print DESC "\n";
 foreach my $key ( #
   sort { $hash->{$b}->{$field} <=> $hash->{$a}->{$field} } #
   keys %{$hash}
 )
 {
  my $value = $hash->{$key};
  #printf( "%s %s\n", $value->{'count'}, $value->{'desc'} );#DEBUGCODE
  print DESC join("\t",$rowCount ,$key, $value->{$field},$value->{'count'}, $value->{'desc'});
  print DESC "\n";
  $rowCount++;
 }
 close(DESC);
} 


=head2 sortAndWriteEnrichment

 Title : sortAndWriteEnrichment
 Usage : sortAndWriteEnrichment( -fileName => 'the output ',
        -hash => 'an hash with the pairs GOID, occurrences',
         -header => 'the header to be written'
     );

 Function: this function starts from an hash with the informations about bp, mf, cc or domains and
   order it from the highest number of occurrences to the lowest. 
   All the pairs sorted are written in an input file. The difference with sortAndWrite is that
   this subroutine is used to order a single argument hash hash->{value}
   Where the counter is incremented using the ++ notations: hash->{value}++

 Returns : a file with the sorted elements

=cut
sub sortAndWriteEnrichment{
 my $fileName=shift;
 my $hash = shift;
 my $header =shift;
 
 my $rowCount = 0;
 
 open(DESC, ">".$fileName);
 print DESC join("\t",@$header);
 print DESC "\n";
 foreach my $key (  #
  sort { $hash->{$b}->{'explev'} <=> $hash->{$a}->{'explev'} }  #
  keys %{$hash}
  )
 {
  my $value = $hash->{$key};
  #printf( "%s %s\n", $value->{'count'}, $value->{'desc'} );
  print DESC join("\t",$rowCount ,$key,$value->{'explev'}, $value->{'desc'});
  print DESC "\n";
  $rowCount++;
 }
 close(DESC);
} 


=head2 sortAndWriteSingle

 Title : sortAndWriteSingle
 Usage : sortAndWriteSingle( -fileName => 'the output ',
        -hash => 'an hash with the pairs GOID, occurrences',
         -header => 'the header to be written'
     );

 Function: sort and write a file with the most enriched GO Terms

 Returns : a file with the sorted elements

=cut
sub sortAndWriteSingle{
 my $fileName=shift;
 my $hash = shift;
 my $header =shift;

 my $rowCount = 0;
 
 open(DESC, ">".$fileName);
 print DESC join("\t",@$header);
 print DESC "\n";
 foreach my $key (  #
  sort { $hash->{$b} <=> $hash->{$a} }  #
  keys %{$hash}
  )
 {
  my $value = $hash->{$key};
  print DESC join("\t",$rowCount ,$key,$value);
  print DESC "\n";
  $rowCount++;
 }
 close(DESC);
} 


=head2 createCDHash

 Title  : createCDHash
 Usage  : createCDHash( -dbDataFolder => 'folder with the databases data',
              -database => 'reference to the proteins database',
              -mySqlUser => 'password of the database',
              -mySqlPass => 'user-id of the database',
              -cdTableHash => 'hash used for the domains info'
       );

 Function: creates the CD Table with each of the domains descriptions. Descriptions
     are taken from a downloaded file. They are 5 columns separated by TAB.
     An hash given in input (cdTableHash) is built using the informations about the domains and
     it will be used to add the informations to the final tables.
 Returns : nothing

=cut
sub createCDHash{
 #Input Variables
 my $dbDataFolder = shift;
 my $database = shift;
 my $mySqlUser = shift;
 my $mySqlPass = shift;
 my ($cdTableHash) =shift;
 
 #CREATING cd HASH
 # PERL DBI CONNECT AND CREATE TABLE
 my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
 if ($dbh->err()) { die "$DBI::errstr\n"; }
 
 print "Creating an hash for the conserved domains....";
 my $select_handle = $dbh->prepare('SELECT * FROM cdTable');
 die "ERROR [$!]: Couldn't prepare queries; aborting" unless defined $select_handle;
 $select_handle->execute();
 
 while (my $res = $select_handle->fetchrow_hashref){
  $$cdTableHash->{$res->{'cdName'}}->{'cd_id'} = $res->{'cd_id'};
  $$cdTableHash->{$res->{'cdName'}}->{'genericId'} = $res->{'genericId'};
  $$cdTableHash->{$res->{'cdName'}}->{'cdName'} = $res->{'cdName'};
  $$cdTableHash->{$res->{'cdName'}}->{'cdShortName'} = $res->{'cdShortName'};
  $$cdTableHash->{$res->{'cdName'}}->{'cdDescription'} = $res->{'cdDescription'};
 }
 
 print "....Completed!\n";
 $dbh->commit();
 $dbh->disconnect();
}


=head2 detectLORFAgreement

 Title : detectLORFAgreement
 Usage : detectLORFAgreement( -orfFile => 'the longest ORF file ',
                              -outFolder => 'path to the outputs folder',    
     );

 Function: The subroutine controls the agreement among strands and frames from results obtained with Swiss-Prot,
            Trembl and Uniref. It reads from a file created with the longest ORF and this kind of informations.
            
            This script behaves differently in this 3 different situations:
            #>comp11465_c0_seq1 -|-|-|DNA2PEP|225|-|2|
            #>comp11465_c0_seq1 Swiss-Prot|Q95V11|-1|2|UniRef90|Q95V11|-1|2|DNA2PEP|225|-|2|
            #>comp11465_c0_seq1 UniRef90|Q95V11|-1|2|DNA2PEP|225|-|2|
            
            In the first case there are no results for Swiss-Prot Trembl (or Uniref)->No evaluation will be done
            In the second both results (SP and Uniref or SP and TR) are present -> the match is double
            In the third case only one is present (SP,TR or Uniref) -> single match
            
            In the case of a single match the indexes are s1 and s2. For a double they are s1 with s3 and s2 with s3.
 Returns : the number of sequences in agreement
=cut
sub detectLORFAgreement{
  my $orfFile = shift;
  my $outFolder = shift;
  
  $orfFile = $outFolder.'/'.$orfFile;
  open(IN,"<$orfFile ");
  
  my $inAgree = 0;
  my $strandORF;
  
  #We initialize with numbers corresponding to the fact that all the elements are present
  my $s1 = 2;
  my $f1 = 3;
  my $s2 = 6;
  my $f2 = 7;
  my $s3 = 10;
  my $f3 = 11;
  
  my $toMatch = 0;
  
  while(my $row = <IN>) {
   #If the line starts with a > then it will save only the string until the first space
   if($row =~ /^\>/){
    $row =~ s/>//g;
    #Separating the name from the data used
    my @pieces = split(/(\Q$separator\E)/,$row);
    my @strandInfos = split('\|',$pieces[2]);
    
      #IF at least one of the results is printed
     if( ($strandInfos[1] ne '-') ){
        if( $strandInfos[0] =~ 'SwissProt' ){#Swiss-prot result is there
          if ( $strandInfos[4] =~ 'DNA2PEP' ){#Only sp is there
            $toMatch = 1;
          }else{#There is something and then DNA2PEP result
            $toMatch = 2;
          }
        }else{
            $toMatch = 1;
        }
              
      if ($toMatch == 1){
        #Change the type of ORF to match - is -1 and + is 1
       if($strandInfos[$s2] eq '-'){
        $strandORF = '-1';
       }else {$strandORF = '1'}
       if( ($strandInfos[$s1] eq $strandORF) and ($strandInfos[$f1] eq $strandInfos[$f2])){
         $inAgree++;
        }
      }elsif ($toMatch == 2){
        #Change the type of ORF to match - is -1 and + is 1
       if($strandInfos[$s3] eq '-'){
        $strandORF = '-1';
       }else {$strandORF = '1'}
        if( ($strandInfos[$s1] eq $strandORF) and ($strandInfos[$s2] eq $strandORF) and
         ($strandInfos[$f1] eq $strandInfos[$f3]) and ($strandInfos[$f2] eq $strandInfos[$f3])){
          $inAgree++;#increase the agreement of sequences
        }
      }#else there is nothing to match
     }
   }
  }
  close(IN);
 return $inAgree;
}

 
 
=head2 extractPercentages

 Title : extractPercentages
 Usage : extractPercentages( -fastaSeqs => 'the sequences file ',
                              -totSeqs => 'number of sequences',
                               -outFile => 'an output file with this statistic'
     );

 Function: extract percentages of ATGC and N and GC content. It uses bioperl to do this. The output file is created in the
  working directory

 Returns : an array of statistics
=cut
sub extractPercentages{
 my $fastaSeqs = shift;
 my $totSeqs = shift; 
 my $outFile = shift;
 
 my $in = Bio::SeqIO->new(-format => 'fasta',
        -file=> $fastaSeqs);
 my ($seqcount,$basecount,$A,$G,$C,$T,$N,$E,$X,$GC);

 open (OUT, ">$outFile") or die "ERROR [$!]: Cannot open $outFile to read... \n";
 my @headerOUT= ("ID\t"."Length\t"."Adenine\t"."Cytosine\t"."Guanine\t"."Thymine\t"."N\t"."GCcontent\t"."Errors\n");
 print OUT @headerOUT;


 my $aPerc = 0 ;
 my $gPerc = 0 ;
 my $tPerc = 0 ;
 my $cPerc = 0 ;
 my $gcPerc = 0;
 my $nPerc = 0 ;
 my $tot_seq_len = 0;
  
 $E=0;
 while ( my $seq = $in->next_seq()) {
  

  my $display= $seq->display_id;
  #my $description= $seq->desc();
  my $str = $seq->seq; # get the sequence as a string
  
  #Erase spaces
  $str =~ s/\s+//g;
  #Count occurrences of A,G,C,T and N
  my $eachA = ($str =~ tr/Aa/Aa/);
  $A += $eachA;
  my $eachC = ($str =~ tr/Cc/Cc/);
  $C += $eachC;
  my $eachG = ($str =~ tr/Gg/Gg/);
  $G += $eachG;
  my $eachT = ($str =~ tr/Tt/Tt/);
  $T += $eachT;
  my $eachN = ($str =~ tr/Nn/Nn/);
  $N += $eachN;
  
  #Counts the occurrences of G and C
  #my $eachGC += ($str =~ tr/GC/gc/);
  #$GC += $eachGC;
  
  #Let's store the overall count
  my $sum = $A + $G + $C + $T + $N; 
  my $Adenine = (($A/$sum)*100);
  my $Cytosine = (($C/$sum)*100);
  my $Guanine = (($G/$sum)*100);

  my $Thymine = (($T/$sum)*100);
  my $NN=(($N/$sum)*100);
  my $basecount = $seq->length;
  
  #The GC content calculation
  my $single_GC = 0;
  $tot_seq_len += $basecount;
  while ($str =~ /[gc]/ig) {
     $GC++;
     $single_GC++;
  }
  #number of occurrences of G close to G (CG) divided per the length of the sequence
  my $GCcontent = ($single_GC/$basecount)*100;
  #RESULTS HAVE BEEN COMPARED WITH endmemo.com!
  
  #print on the file the information for the current sequence
  $E = ($basecount) - ($A + $G + $C + $T + $N);
  print OUT $seq->id."\t".$basecount."\t".$Adenine."\t".$Cytosine."\t".$Guanine."\t".$Thymine."\t".$NN."\t".$GCcontent."\t".$E."\n";

  $aPerc = $aPerc + $Adenine;
  $gPerc = $gPerc + $Guanine;
  $tPerc = $tPerc + $Thymine;
  $cPerc = $cPerc + $Cytosine;
  $gcPerc = $gcPerc + $GCcontent;
  $nPerc = $nPerc + $NN;
  
  $A=0;
  $G=0;
  $C=0;
  $T=0;
  $N=0;
  $E=0;

}
 #overall GC content
 my $tot_GC_cont = ($GC / $tot_seq_len) * 100;
 
 
 return ($aPerc/$totSeqs,$gPerc/$totSeqs,$tPerc/$totSeqs,$cPerc/$totSeqs,$nPerc/$totSeqs,$tot_GC_cont); 
}


1;
