#############################programs_out_parser.pm#################
# This script uses the information taken by the programs (BLAST, lncRNA predicttion, DNA2PEP) to
# create some GFF3 files and a database to access those informations. At the same time an hash that contain them all is created
# to access rapidly the data when a final table is written.
############################################################

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


package LIB::programs_out_parser;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw(execute_programs_out_parser);


#GENERAL UTILITY MODULES
use strict;
use warnings;
use Data::Dumper;#To print hash dumper
use Getopt::Long;#To control the input
use Storable;
use Time::HiRes qw( time );#To keep in count the times of Annocript
use Error qw(:try);#To alert about errors
use IO::Handle;#To immediately print with autoflush 

#SEQUENCE MANAGEMENT MODULES
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;#To use gff format
use Bio::SearchIO;#To use BioPerl sequence management
use Bio::SeqFeature::Generic;#To access sequences features
use Bio::Root::Exception;

#SQL DATA MANAGEMENT MODULES
use DBI;

#FILE MANAGEMENT MODULES
use LWP::UserAgent; #Ftp management
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);#Compression and uncompression of data


#Annocript libraries
use LIB::output_and_stats	qw( lncRNA_coding_fasta_file_generator lncRNA_fasta_file_generator	
															extractAnnotated extractStatistics );

use LIB::files_management qw( extract_name checkOutputsPresence delete_file);

use LIB::mysql_db_management	qw( drop_database destroy_db );

#Using a library for standard utilities								
use LIB::std_lib qw(print_array);	

use LIB::programs_management qw(try_exec_command);	


my $configHash;# Main hash with configuration variables
my $gffFolder;#Path to the folder with gff files
my $outFolder;#Path to the folder with outputs
my $dataFolder;#Path to the data folder of the session
my $sessionFolder;#Path to the session folder
my $dbDataFolder;#Path to the folder with the db data
my $logFolder;#Path to the folder with the log files

my $database;#database used
my $platform;#platform of database used (mysql)
my $host;#host to connect the database
my $port;#Port for db connection
my $mySqlUser;#user id for db connection
my $mySqlPass; #password for db connection
my $dsn;#dsn for db connection
my $dsnGff;#dsn for gff db connection
my $gffDB;#name of the GFF db
my $dbargs;#all the args of db connection in a variable

	 
my $outHash;#An hash to contain all the element to be written in the output file	
 
my $separator;#separator of elements in the final table
my $os_tax;#organism or taxonomy: a variable to define what is in the database
my $os_taxTable;#organism or taxonomy: a variable to define which table is in the database
my $view;#A view to access rapidly data in the database

my @tr_ids = ();#Array which contains all the transcript names to analyze

my @trWords = ('uncharacterized','predicted protein','hypothetical','putative uncharacterized');#Words to match to be a non coding

#HEADERS 
#These are the headers of the output table. Here you can add fields if you need
my @transcriptsCharHeader = qw(TranscriptName
                            TransLength);
my @blastxSPTRHeader  = qw(HSPNameSP
                            HSPLengthSP
                            HSPEvalueSP
                            HITLengthSP
                            QCoverageSP
                            HCoverageSP
                            StrandSP
                            FrameSP
                            DescriptionSP
                            OSNameSP
                            EnzymeIds
                            EnzymeDescs
                            PwLev1
                            PwLev2
                            PwLev3
                            HSPNameTR
                            HSPLengthTR
                            HSPEvalueTR
                            HITLengthTR
                            QCoverageTR
                            HCoverageTR
                            StrandTR
                            FrameTR
                            DescriptionTR 
                            OSNameTR);
my @blastxSPUfHeader = qw(HSPNameSP
                            HSPLengthSP
                            HSPEvalueSP
                            HITLengthSP
                            QCoverageSP
                            HCoverageSP
                            StrandSP
                            FrameSP
                            DescriptionSP
                            OSNameSP
                            EnzymeIds
                            EnzymeDescs
                            PwLev1
                            PwLev2
                            PwLev3
                            HSPNameUf
                            HSPLengthUf
                            HSPEvalueUf
                            HITLengthUf
                            QCoverageUf
                            HCoverageUf
                            StrandUf
                            FrameUf
                            DescriptionUf
                            Taxonomy);
my @protGOHeader = qw(BPId BPDesc MFId MFDesc CCId CCDesc);  
my @domGOHeader = qw(domBPId domBPDesc domMFId domMFDesc domCCId domCCDesc);  
my @rpstblastnHeader = qw(CDName
                            CDStartEnd
                            CDEvalue
                            CDDesc);
my @blastnHeader = qw(OtherNCName
                            OtherNCEvalue
                            OtherNCDesc);
my @orfHeader = qw( LongOrfLength
                            LongOrfStrand
                            LongOrfFrame);
my @ncHeader = qw(ProbToBeNonCoding);
my @heuristicHeader = qw(lncRNA4Annocript);
my @seqHeader = qw(Sequence);


#IMPORTANT: BEFORE TO ADD LOOK AT THE DOCUMENTATION OF THE SOFTWARE!!!!
            
my $finalTabHeader;#This will be the chose one. It will be valued once that the db has been decided
my $ridFinalTabHeader;

  
=head2 execute_programs_out_parser

 Title : execute_programs_out_parser
 Usage : execute_programs_out_parser( -configHash => 'configHash',
					);

 Function: Takes in input the configuration hash and calls all the routines needed to print the output and the gff files.
						Finally it prints the statistics.

 Returns : some output and statistics files

=cut	 
sub execute_programs_out_parser{
  
	#This is an hash table with all the elements read in the config.txt file
	$configHash = shift;
  
  STDOUT->autoflush(1);#This makes STDOUT hot in the sense that everything will be print immediately

  #Set some global variables by using the configuration hash
	$sessionFolder = $configHash->{'sessionFolder'};
	$gffFolder = $sessionFolder."/".$configHash->{'gffFolder'};
	$dataFolder = $sessionFolder."/".$configHash->{'dataFolder'};
	$outFolder = $sessionFolder."/".$configHash->{'outFolder'};
	$dbDataFolder = $sessionFolder."/".$configHash->{'dataFolder'};
	$logFolder = $sessionFolder."/".$configHash->{'logFolder'};
	
	#Database parameter configuration
	$database = $configHash->{'database'};
	$mySqlUser = $configHash->{'mySqlUser'};
	$mySqlPass = $configHash->{'mySqlPass'};
	$platform = $configHash->{'platform'};
	$host = $configHash->{'host'};## MySQL Host
	$port = $configHash->{'port'};
	$dbargs = {AutoCommit => 0, PrintError => 1};
	$dsn = "dbi:$platform:$database:$host:$port";	
	$gffDB = $configHash->{'gffDB'}; 
  $dsnGff = "dbi:$platform:$gffDB:$host:$port;mysql_local_infile=1";

  #Set the separator
	$separator = $configHash->{'separator'};
	
  #Taking the time to count
	my $globalStart = time;
	my $partTime = time;
	my $partDuration = undef;

  #First an existing database with that name will be erased
  if ($configHash->{'destroy_gffdb'} eq "YES"){
	  drop_database($configHash->{'gffDB'},$platform,$host,$port,$mySqlUser,$mySqlPass) or die "Unable to destroy ".$configHash->{'gffDB'}.": ?";
	}
	
	#Storing the fasta sequences in the GFF database
	fasta_2_gffdb($configHash->{'fastaSeqs'},$configHash->{'gffDB'},"Y");
	@tr_ids = collectFastaNames($configHash->{'fastaSeqs'});#fasta_2_array_of_ids($configHash->{'fastaSeqs'});
	
	$partDuration = time-$partTime;
	print "Finished storing sequences in ".$partDuration." seconds\n";
	$partTime = time;
		
 
  if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){
    $os_tax='os';
    $os_taxTable = 'OSType';
    $view = 'uniprotkbViewSeq';
    print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE BLASTX (BLASTP) OUTPUT WITH SWISS PROT<<<<<<<<<<<<<<\n";#LOG
    #Extracting information from the BLAST out for swiss-prot and creating GFF3 file						
    if ( multiBlastout_2_gff3AndHash($configHash->{'blastxSPOut'},$configHash->{'blastxSPOutType'},$configHash->{'blastxSPGFF3FileName'},"blastx",1,$configHash->{'database'},'sp','OSType','os','uniprotkbViewSeq') > 0){
      $configHash->{'blastxSPOutPresent'} = 'YES';#A variable to know if this output has been created
      $partDuration = time - $partTime;
      print "Finished storing output from blastx against SwissProt in\t$partDuration\n";
      $partTime = time;
     
      if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'blastxSPGFF3FileName'}) ){
        gff3_2_db($configHash->{'blastxSPGFF3FileName'},$configHash->{'gffDB'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
        $partDuration = time - $partTime;
        print "Finished generation of GFF3 file in $partDuration\n";
        $partTime = time;
      }
      
    }else{$configHash->{'blastxSPOutPresent'} = 'NO';}
    
    print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE BLASTX (BLASTP) OUTPUT WITH TREMBL<<<<<<<<<<<<<<\n";#LOG
    #Extracting information from the BLAST out for trembl and creating GFF3 file					
    if ( multiBlastout_2_gff3AndHash($configHash->{'blastxTROut'},$configHash->{'blastxTROutType'},$configHash->{'blastxTRGFF3FileName'},"blastx",1,$configHash->{'database'},'tr','OSType','os','uniprotkbViewSeq') > 0 ){
      $configHash->{'blastxTROutPresent'} = 'YES';#A variable to know if this output has been created
      $partDuration = time - $partTime;
      print "Finished storing output from blastx against TrEMBL in\t$partDuration\n";
      $partTime = time;
      
      if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'blastxTRGFF3FileName'}) ){
        $partDuration = time - $partTime;
        gff3_2_db($configHash->{'blastxTRGFF3FileName'},$configHash->{'gffDB'});
        print "Finished generation of GFF3 file in $partDuration\n";
        $partTime = time;  
      }
      
    }else{$configHash->{'blastxTROutPresent'} = 'NO';}
 }
 
 if ( $configHash->{'dbInUse'} eq 'uniref' ){
   $os_tax='os';
    $os_taxTable = 'OSType';
    $view = 'unirefViewSeq';
    
    print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE BLASTX (BLASTP) OUTPUT WITH SWISS PROT<<<<<<<<<<<<<<\n";#LOG
    #Extracting information from the BLAST out for swiss-prot and creating GFF3 file						
    if ( multiBlastout_2_gff3AndHash($configHash->{'blastxSPOut'},$configHash->{'blastxSPOutType'},
    $configHash->{'blastxSPGFF3FileName'},"blastx",1,$configHash->{'database'},'sp','OSType','os','uniprotkbViewSeq') > 0){
      $configHash->{'blastxSPOutPresent'} = 'YES';#A variable to know if this output has been created
      $partDuration = time - $partTime;
      print "Finished storing output from blastx against SwissProt in\t$partDuration\n";
      $partTime = time;
     
      if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'blastxSPGFF3FileName'}) ){
        gff3_2_db($configHash->{'blastxSPGFF3FileName'},$configHash->{'gffDB'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
        $partDuration = time - $partTime;
        print "Finished generation of GFF3 file in $partDuration\n";
        $partTime = time;
      }
      
    }else{$configHash->{'blastxSPOutPresent'} = 'NO';} 
    
   $os_tax='tax';
   $os_taxTable = 'taxType';
   print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE BLASTX (BLASTP) OUTPUT WITH UNIREF<<<<<<<<<<<<<<\n";#LOG
   #Extracting information from the BLAST out for swiss-prot and creating GFF3 file						
    if ( multiBlastout_2_gff3AndHash($configHash->{'blastxUnirefOut'},$configHash->{'blastxUnirefOutType'},
    $configHash->{'blastxUnirefGFF3FileName'},"blastx",1,$configHash->{'database'},'uniref','taxType','tax','unirefViewSeq') > 0){
      $configHash->{'blastxUnirefOutPresent'} = 'YES';#A variable to know if this output has been created
      $partDuration = time - $partTime;
      print "Finished storing output from blastx against UniRef in\t$partDuration\n";
      $partTime = time;
     
     if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'blastxUnirefGFF3FileName'}) ){
       gff3_2_db($configHash->{'blastxUnirefGFF3FileName'},$configHash->{'gffDB'},$configHash->{'mySqlUser'},$configHash->{'mySqlPass'});
       $partDuration = time - $partTime;
       print "Finished generation of GFF3 file in $partDuration\n";
       $partTime = time;
      }
      
    }else{$configHash->{'blastxUnirefOutPresent'} = 'NO';}
 }
	#print "all has been done in ".$partDuration." seconds\n";
	#$partTime = time;

  print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE RPS BLAST OUTPUT<<<<<<<<<<<<<<\n";#LOG
	#converting a rpstblastn output to a gff3 file
	if ( multiBlastout_2_gff3AndHash($configHash->{'rpstblastnOut'},$configHash->{'rpstblastnOutType'},
    #$configHash->{'rpstblastnGFF3FileName'},"rpstblastn",2,undef,'none',undef,undef,'cdViewGO') > 0){
    $configHash->{'rpstblastnGFF3FileName'},"rpstblastn",2,$configHash->{'database'},'none',undef,undef,'cdViewGO','printGFF') > 0){

    $configHash->{'rpstblastnOutPresent'} = 'YES';#A variable to know if this output has been created
    $partDuration = time - $partTime;
    print "Finished storing output from rpsblast against domains in\t$partDuration seconds\n";
    $partTime = time;
    
    $partDuration = time - $partTime;
    gff3_2_db($configHash->{'rpstblastnGFF3FileName'},$configHash->{'gffDB'});
    print "Finished generation of GFF3 file in $partDuration seconds\n";
    $partTime = time;
  }else{$configHash->{'rpstblastnOutPresent'} = 'NO';}
  
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE BLASTN (TBLASTN) OUTPUT<<<<<<<<<<<<<<\n";#LOG
	#converting a blastx output to a gff3 file
	if ( multiBlastout_2_gff3AndHash($configHash->{'blastnOut'},$configHash->{'blastnOutType'},
    $configHash->{'blastnGFF3FileName'},"blastn",2,undef,'none',undef,undef,undef) > 0){
    $configHash->{'blastnOutPresent'} = 'YES';#A variable to know if this output has been created
    $partDuration = time - $partTime;
    print "Finished storing output from blastn against ncRNAs in\t$partDuration seconds\n";
    $partTime = time;
    
    if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'blastnGFF3FileName'}) ){
      #Storing the gff file in the GFF database
      gff3_2_db($configHash->{'blastnGFF3FileName'},$configHash->{'gffDB'});
      $partDuration = time - $partTime;
      print "Finished generation of GFF3 file in ".$partDuration." seconds\n";
      $partTime = time;
    }
   
  }else{ $configHash->{'blastnOutPresent'} = 'NO';}
  
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE lncRNA prediction program OUTPUT<<<<<<<<<<<<<<\n";#LOG
	#NON CODING 
	#if ( portraitout_2_gff3AndHash($configHash->{'lncRNAPredOut'},$configHash->{'nonCodingAlgorithm'},$configHash->{'nonCodingGFF3FileName'}) > 0){
	if ( lncRNAPredOut_2_gff3AndHash($configHash->{'lncRNAPredOut'},$configHash->{'nonCodingAlgorithm'},$configHash->{'nonCodingGFF3FileName'}) > 0){
     $configHash->{'lncRNAPredOutPresent'} = 'YES';#A variable to know if this output has been created
    $partDuration = time - $partTime;
    print "Finished storing output from CPC2 in $partDuration seconds\n";
    
    if ( ($configHash->{'useGFFDB'} eq 'YES')  and (-e $configHash->{'nonCodingGFF3FileName'}) ){
      $partTime = time;
      gff3_2_db($configHash->{'nonCodingGFF3FileName'},$configHash->{'gffDB'});
      $partDuration = time - $partTime;
      print "Finished generation of GFF3 file in $partDuration seconds\n";
      $partTime = time;
    }
  }else{$configHash->{'lncRNAPredOutPresent'} = 'NO';}
  
  print "\n(".scalar(localtime).") >>>>>>>>>>>>>EXTRACTION OF INFO FROM THE DNA2PEP OUTPUT<<<<<<<<<<<<<<\n";#LOG
	#LONGEST ORF
	if ( dna2pepout_2_gff3AndHash($configHash->{'dna2pepOut'},$configHash->{'orfAlgorithm'},$configHash->{'orfGFF3FileName'}) > 0){
    $partDuration = time - $partTime;
    print "Finished storing output from dna2pep in\t".$partDuration." seconds\n";
    $configHash->{'dna2pepOutPresent'} = 'YES';#A variable to know if this output has been created
    if ( ($configHash->{'useGFFDB'} eq 'YES') and (-e $configHash->{'orfGFF3FileName'}) ){
      $partTime = time;
      gff3_2_db($configHash->{'orfGFF3FileName'},$configHash->{'gffDB'});
      $partDuration = time - $partTime;
      print "Finished generation of GFF3 file in ".$partDuration." seconds\n";
      $partTime = time;
    }
  }else{$configHash->{'dna2pepOutPresent'} = 'NO';}
  
  
  #Save the hash just created in a file to use it later in the print  
  my $outHashFile = $dataFolder."/".$configHash->{'outHashFile'};
	saveOutHash($outHashFile); # Print the dumper on a file
  
  #Here we store the complete hash in a file
  #This operation is done to have the hash saved and in a future Annocript can  
  #also separate the tasks
  #my $outHashFile = $dataFolder."/".$configHash->{'outHashFile'};
  #store $outHash, $outHashFile;
  
  #printSequencesID("gffDB","mysql_dev","dEvEl0pEr");#DEBUGCODE
  
  #Print all the info in the final table
	print "\n(".scalar(localtime).") >>>>>>>>>>>>>GENERATION OF THE TABULAR OUTPUT<<<<<<<<<<<<<<\n";#LOG
  printAllInfo($configHash->{'gffDB'},$configHash->{'outFileName'},$configHash->{'outFiltered'},$outHashFile);#$configHash->{'storeBlastx'},$outHashFile);
	$partDuration = time - $partTime;
	print "Finished printing the complete annotation in ".$partDuration." seconds\n";
	$partTime = time;
  
  #Check if the dna2pep and portrait outputs are present. If they are, then some other file will be created
  $configHash->{'dna2pepOutPresent'} = checkOutputsPresence($configHash->{'dna2pepOut'},$outFolder);

  
  #Print a FASTA file with all the aminoacidic sequences
  if ($configHash->{'dna2pepOutPresent'} eq 'YES' ){
    #creation of the ORF file
    printORFFastaFile( $configHash->{'ORFFastaFileName'},$outFolder);
    $partDuration = time - $partTime;
    print "Finished writing orf informations in ".$partDuration." seconds\n";
    $partTime = time;
  }

   #If user chose to not create the gff file its use can finish here. 
  #No access will be provided to the GFF Database
  if ( $configHash->{'useGFFDB'} eq 'NO' ){
    print "You don't need any GFF database so ".$configHash->{'gffDB'}." will be destroyed.\n";
    drop_database($configHash->{'gffDB'},$platform,$host,$port,$mySqlUser,$mySqlPass);
  }
 
  #Calculate the time needed for all the process
	my $globalDuration = time - $globalStart;
  print "\n".scalar(localtime)." Annocript has finished the creation of the output table in ".$globalDuration." seconds\n";
}



=head2 dna2pepout_2_gff3AndHash

 Title : dna2pepout_2_gff3AndHash
 Usage : dna2pepout_2_gff3AndHash( -file => 'dna2pepout_2out',
					);

 Function: creates a gff3 file from a dna2pep output and an Hash containing all the data. 

 Returns : a gff3 file

=cut
sub dna2pepout_2_gff3AndHash{
	my $dna2pepOutFile = shift;
	my $algorithmType = shift;
	my $gff3FileName = shift;

	my $return = 0;#If it's 1 then the subroutine has been executed
  #If the output from dna2pep is present then extract the output
  
  $dna2pepOutFile = $outFolder."/".$dna2pepOutFile;
  if ( -e $dna2pepOutFile ){
		#Open the output from the program
		open (DNA2PEP_OUT,"<$dna2pepOutFile") || print ( "\nWARNING: Cannot open dna2pep out file for read!!!" );
		
		#If the user wants the GFF output
		if ($configHash->{'printGFFOutput'} eq 'YES'){
			$gff3FileName = $gffFolder."/".$gff3FileName;
			open (GFF_FILE,">$gff3FileName") || die ( "\nERROR [$?]: Cannot open gff file to write: ?" );
			print GFF_FILE "##gff-version 3 \n";
		}
		#Inizialization of variables
		my $orfLength = undef;
		my $featId = 1;
		my $strand = undef;
		my $frame = undef;
		
    while(<DNA2PEP_OUT>){
      $return = 1;
      my $data = $_;
      chomp $data;
        
      # For example we show here a string with all the result of dna2pep
      # comp177849_c0_seq1_rframe-3_ORF	FSQLFPSTILPKLY	.......m..........m /orf_mode="none"; /dna="AGCCTAAGTTCC"
      # First we split all the line based on the tab
      my @fields = split("\t", $data);
      
      #extract the frame that is written after the rframe tag
      $fields[0] =~ s/_rframe(.+)_ORF//g;
      my $rFrame = $1;
      
      #now if it is positive we put a + otherwise a -
      if( $rFrame < 0 ){
        $strand = "-";
      }else{$strand = "+";}
      
      #and, in any case, take the abs value for the frame	
      #DNA2PEP prints the frame as 1,2 or 3 while blast uses 0,1,2
      #we here reduce the number of dna2pep for uniformity in all the output files	
      $frame = abs($rFrame)-1;
        
      my $id = $fields[0];
      if (grep {$id} @tr_ids){
				#$id =~ s/[^A-Za-z0-9\-\(\)\_\-\.\/]/_/g;
				###
				#The orf and its length is computed taking the first element of the first separation
				my $orf = $fields[1];
				my $orfLength = length($orf); 
				$orf =~ s/[^A-Z]//g;#CHECK THIS IF THE FUNCTION WORKS WITHOUT
				
						#If the user wants the GFF output
				if ($configHash->{'printGFFOutput'} eq 'YES'){	
					#print GFF_FILE "##sequence-region ".$id." 1 ".$orfLength."\n";
					print GFF_FILE $id."\t".$algorithmType."\tORF\t";
					print GFF_FILE "10\t20\t0\t".$strand."\t".$frame."\t";
					print GFF_FILE "l=".$orfLength.";os=".$strand.";of=".$frame.";orf=".$orf."\n";#"ID=".$featId.";
				}
				$featId++;
					
				#Adding Value to global Hash
				$outHash->{$id}->{'dna2pep'}->{'l'} = $orfLength; 
				$outHash->{$id}->{'dna2pep'}->{'os'} = $strand; 
				$outHash->{$id}->{'dna2pep'}->{'of'} = $frame;
				$outHash->{$id}->{'dna2pep'}->{'orf'} = $orf;
		  }  
            
      #print "id: ".$id." - frame: ".$frame." - strand: ".$strand." orf: ".$orf." - len: ".$queryLength."\n";#DEBUGCODE
    }
    close(DNA2PEP_OUT);
  }else{ print("Output from dna2pep has not been extracted because $dna2pepOutFile has not been found\n");}
  
  #If the user wants the GFF output
	if ($configHash->{'printGFFOutput'} eq 'YES'){
		close(GFF_FILE);
		print "Your file $gff3FileName has been correctly written!\n ";
	}
  return $return;
}

=head2 lncRNAPredOut_2_gff3AndHash

 Title : lncRNAPredOut_2_gff3AndHash
 Usage : lncRNAPredOut_2_gff3AndHash( -file => 'dna2pepout_2out',
					);

 Function: creates a gff3 file from the output of long non coding prediction and an hash with same informations. 
          Example of CPC2 output: >comp165_c0_seq1 :1: 0.948593 0.0514073
  #ID	           transcript_length	peptide_length	Fickett_score	pI	       ORF_integrity	coding_probability	label
	comp53_c1_seq1 582	              43	             0.34424	  5.95758056641	    1	             0.0131439	  noncoding

 Returns : a gff3 file

=cut
sub lncRNAPredOut_2_gff3AndHash{
	my $lncRNAPredOutFile = shift;
	my $algorithmType = shift;
	my $gff3FileName = shift;
	 
	$lncRNAPredOutFile = $outFolder."/".$lncRNAPredOutFile ;
  my $return = 0;#If it's 1 then the subroutine has been executed
	my $id_ind = 0;
	my $cod_pot_ind = 6;
	my $peplen_ind = 2;
  if ( -e $lncRNAPredOutFile ){
		
		open (IN_FILE,"<$lncRNAPredOutFile") or die ( "\nERROR [$?]: Cannot open $lncRNAPredOutFile: ?" );
		 #If the user wants the GFF output
		if ($configHash->{'printGFFOutput'} eq 'YES'){
			$gff3FileName = $gffFolder."/".$gff3FileName;
			open (GFF_FILE,">$gff3FileName") or die ( "\nERROR [$?]: Cannot open $gff3FileName to write: ?" );
			print GFF_FILE "##gff-version 3 \n";
	  }
    my $featId = 1;
  

    while(<IN_FILE>){
      $return = 1;
      my $data = $_;
      chomp $data;
      
      #Split the name from the information. This is the only clean way (separate with :1:)
      my @pieces = split('\t',$data);
      #print "name:$pieces[0], the rest: $pieces[1]\n";#DEBUG CODE
      
     #The first element is the id.
      my $name = $pieces[$id_ind];
      
      #print "the name is: $name\n";#DEBUG CODE
      
      if (grep {$name} @tr_ids){
				#Get the coding potential
				my $coding_potential = $pieces[$cod_pot_ind];
				my $peplen = $pieces[$peplen_ind];
				#print "$coding_potential - $peplen\n";#DEBUG CODE
				
				 #If the user wants the GFF output
				if ($configHash->{'printGFFOutput'} eq 'YES'){
					 #Write the GFF3 file 
					 print GFF_FILE $name."\t".$algorithmType."\t"."ncRNA\t";#WE CHOOSE ncRNA since they are sequences evaluated for this
					 print GFF_FILE "10\t20\t$coding_potential\t.\t.\t";#10-20 ARE INVENTED!!!!
					 print GFF_FILE "ID=".$featId.";codpot=$coding_potential;peplen=$peplen\n";#codpot=percentage of be coding - peplen=peptide length
				}
				 $featId++;	
				 
				 #Adding Value to global hash
				 $outHash->{$name}->{'codpot'} = $coding_potential;
			}
     }
    close(IN_FILE);
    #If the user wants the GFF output
		if ($configHash->{'printGFFOutput'} eq 'YES'){
			close(GFF_FILE);
			print("GFF FILE for coding/non-coding has been created!\n");
		}
  }else{ print("GFF FILE for coding/non-coding has not been created. $lncRNAPredOutFile is not present.\n");}
  
 return $return;
}

#=head2 portraitout_2_gff3AndHash

 #Title : portraitout_2_gff3AndHash
 #Usage : portraitout_2_gff3AndHash( -file => 'dna2pepout_2out',
					#);

 #Function: creates a gff3 file from the portrait output and an hash with same informations. 
          #Example of Portrait output: >comp165_c0_seq1 :1: 0.948593 0.0514073

 #Returns : a gff3 file

#=cut
#sub portraitout_2_gff3AndHash{
	#my $portraitOutFile = shift;
	#my $algorithmType = shift;
	#my $gff3FileName = shift;
	 
	#$portraitOutFile = $outFolder."/".$portraitOutFile ;
  #my $return = 0;#If it's 1 then the subroutine has been executed
	
  #if ( -e $portraitOutFile ){
		
		#open (PORTRAIT_OUT,"<$portraitOutFile") or die ( "\nERROR [$?]: Cannot open $portraitOutFile: ?" );
		 ##If the user wants the GFF output
		#if ($configHash->{'printGFFOutput'} eq 'YES'){
			#$gff3FileName = $gffFolder."/".$gff3FileName;
			#open (GFF_FILE,">$gff3FileName") or die ( "\nERROR [$?]: Cannot open $gff3FileName to write: ?" );
			#print GFF_FILE "##gff-version 3 \n";
	  #}
    #my $queryLength = 10;#INVENTED NUMBER!!!!
    #my $featId = 1;
  

    #while(<PORTRAIT_OUT>){
      #$return = 1;
      #my $data = $_;
      #chomp $data;
      
      ##Split the name from the information. This is the only clean way (separate with :1:)
      #my @pieces = split(':\d:',$data);
      ##print "name:$pieces[0], the rest: $pieces[1]\n";#DEBUG CODE
      
      ##Take only the first piece of the name
      #my @splitName = split(" ",$pieces[0]);
      #my $name = $splitName[0];
      ##The first element is the id. We erase the > in front of it
      #$name =~ s/\>//;
      ##print "And after the split name is: $name\n";#DEBUG CODE
      
      #if (grep {$name} @tr_ids){
				##Separate the coding value from the non coding one
				#my @splitNCInfo = split(" ",$pieces[1]);
				##print "$splitNCInfo[0] - $splitNCInfo[1]\n";#DEBUG CODE
				
				 ##If the user wants the GFF output
				#if ($configHash->{'printGFFOutput'} eq 'YES'){
					 ##Write the GFF3 file 
					 #print GFF_FILE $name."\t".$algorithmType."\t"."ncRNA\t";#WE CHOOSE ncRNA since they are sequences evaluated for this
					 #print GFF_FILE "10\t20\t".$splitNCInfo[1]."\t.\t.\t";#10-20 ARE INVENTED!!!!
					 #print GFF_FILE "ID=".$featId.";pc=".$splitNCInfo[0].";pnc=".$splitNCInfo[1]."\n";#pc=percentage of be coding - pnc=perc to be non coding
				#}
				 #$featId++;	
				 
				 ##Adding Value to global hash
				 #$outHash->{$name}->{'portrait'} = $splitNCInfo[1];
			#}
     #}
    #close(PORTRAIT_OUT);
    ##If the user wants the GFF output
		#if ($configHash->{'printGFFOutput'} eq 'YES'){
			#close(GFF_FILE);
			#print("GFF FILE for coding/non-coding has been created!\n");
		#}
  #}else{ print("GFF FILE for coding/non-coding has not been created. $portraitOutFile is not present.\n");}
  
 #return $return;
#}


=head2 multiBlastout_2_gff3AndHash
 Title : multiBlastout_2_gff3AndHash
 Usage : multiBlastout_2_gff3AndHash( -file => 'blastout_hspformatted',
					);

 Function: creates a gff3 file and can use different blastx outputs taken from a folder. 
					 If the BLAST program is ran by Annocript, it'll be produced only one out file. Else, if the user
					wants to put its output files in a folder he can do it by changing the names with blastxOut1, blastxOut2,..etc.
					It extracts informations about the hsp and writes on a uniqe file GFF3 according to the standard
          reference format.
          The function works for each of the three blasts: BLASTX (BLASTP), BLASTN (TBLASTN), RPSTBLASTN (RPSBLAST)
          
          To extract informations from different files it needs the name of the output. Then it will add
          iteratively a successive number and searches for that file. If it is found, than it is analyzed.
          
          The database proteins is accessed to extract necessary informations but only for BLASTX they are present
          so I decided to give the variable for the database in input and only when it is present the function access to it.

 Returns : a gff3 file

=cut
sub multiBlastout_2_gff3AndHash{
	
	my $fileName = shift;#Name of the file in input
	my $formatFile = shift;#Format (BLAST, BLASTXML for this one, actually it does not work)
	my $gff3FileName = shift; 
	my $algorithmType = shift;	
	my $indexID = shift;#This is the index, in the name given by the blast output, where can be found the Id
	my $database = shift;#If used this will be the name of the database
	my $uniprotDBType = shift;#sp, tr or uniref
  my $os_taxTable = shift;
  my $os_tax = shift;
  my $view = shift;#The db view to access. It is uniref or uniprotkb based
  my $overridePrintGFF = shift;#used to override the configuration file. Used only for RPS blast
  
  my $return = 0;#If it's 1 then the subroutine has been executed
  
  #Changing the behaviour given in the configuration file when is present the variable
  #overridePrintGFF. Used only for rpsblast output which will be used to create the db GFF
	if (defined  $overridePrintGFF){
		if ($overridePrintGFF eq 'printGFF' )	{
				$configHash->{'printGFFOutput'} = 'YES';
		}
	}
	#Opening database if needed 
	my $dbh;
	if (defined $database){
		# PERL DBI CONNECT AND CREATE TABLE
		$dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1, AutoCommit => 0 });
	}
	
	if ($configHash->{'printGFFOutput'} eq 'YES'){
		$gff3FileName = $gffFolder."/".$gff3FileName;
		#Writing the header of the GFF file only once 
		open (GFF_FILE,">$gff3FileName") || die ( "\nERROR [$?]: Cannot open gff file to write: ?" );
		print GFF_FILE "##gff-version 3 \n";
	}
	
	my $hspId;#Initialize the index of the HSP
  my $hitId;#Initialize the index of the HIT
	my $queryCounter=1;#Inizialization of the query counter
	
	#Here we create a log file where to write informations on 
  #what are the missing information in the headers of TREMBL and SWISSPROT
  my $logFile = $logFolder."/".$uniprotDBType."_".$configHash->{'genesMissingInDB'};
  print "A log of the operation will be written in $logFile\n";
  open (LOG, ">$logFile") 
  or print "ERROR: Cannot create $logFile for missing genes in database. Maybe you need more permissions.\n";
	
	#This is a flag used to understand if there are discordant data in the output with the database
	my $discordantData = 0; 
								
	#It checks	in the output folder for the folder corresponding to the
	#current BLAST then opens iteratively all the files inside this folder.
	my $outputsDir = $outFolder."/".$algorithmType;
	if ( -d $outputsDir){
		#print "DEBUGCODE: Entered since the folder exists..\n";
		opendir DIR, $outputsDir or die "ERROR [$?]: cannot open dir $outputsDir: ? \n";
		my @blastOuts= readdir DIR;
		closedir DIR;
		
		my $numOutFiles = scalar(@blastOuts)-2;#This is necessary to not count . and .. folders
		if ($numOutFiles >0 ){
    print "\nYou have $numOutFiles output files inside the folder $outputsDir...\n ";
		print "Blast output extraction - Algorithm: ".$algorithmType." - Query: ";
			#print "DEBUGCODE: Entered since the outputs are there\n";
			my $index;
			for ($index = 1; $index <= $numOutFiles; $index++){
				my $blastOut = $outputsDir."/".$fileName.$index;
				 #print "DEBUGCODE: trying to access $blastOut\n";
				#Here we access the blast out only if it corresponds to the specific algorithm
				if ( (-e $blastOut) && ($blastOut =~ /$fileName/)){
          $return = 1;#This variable is set to 1 when at least one result is found
					print " \n..Reading ($formatFile) from: $blastOut..\n ";
					my $in = new Bio::SearchIO(-format => $formatFile, -file => $blastOut );
         
					#iterate on results...
					try{
						while( my $result = $in->next_result ) {
							             
              #The Name position change depending from the type of BLAST
							my $qName = $result->query_name;#"";
							#The query name should be present inside the array of transcripts names, so that, if you are
							#using a reduced fasta file, not all the output will be read
							if (grep {$qName} @tr_ids){
																
								$hitId = 1;#Initialize the index of the HIT
								#Print each 1 thousand sequences
								if ( $queryCounter % 1000 == 0 ){
									print "-> $queryCounter ";
								}
								$queryCounter++;
	 			

								#iterate on hits
								while( my $hit = $result->next_hit ) {

									$hspId = 1;#Initialize the index of the HSP
									#The accession sometime is simply the ID but other times it can be separated by pipe or semicolon
									#     293928|20291|agaripopsis
									# BLASTN
									#   gi|225637499|ref|NR_003287.2| Homo sapiens RNA, 28S ribosomal  (FOR RFAM)
									#    RF01059;mir-598;CABZ01027388.1/60-155 (A blastn accession)
									# Here we take only the first element
									my $seqId = $hit->accession;
									my @args;
									if ($hit->accession =~ /\|/){
										@args = split (/\|/, $hit->accession);
										$seqId = $args[1];
													}
									if ($hit->accession =~ /;/){
										@args = split (";", $hit->accession);
										$seqId = $args[2];
									}
									if ( $algorithmType eq "rpstblastn"){
										$seqId = (split(",",$hit->description))[0];
									}
									 
									#Here we open an access to the database to collect informations belonging to the HSPs
									#print "seqID: $seqId\n";
									my @res = ();
									my $select_handle;
									#Extracts information from UNIPROT only if $database is defined (i.e. blastx)	and only if it is the first hit. 
									if ( (defined $database) and ($algorithmType eq "blastx")){ 
										#If the user wants to print the GFF files the query should be performed for each HIT
										#otherwise it is executed only for the first hit
										 if ( ($configHash->{'printGFFOutput'} eq 'YES') or ($hitId == 1) ) {
												#fetch informations from the database
												my $query = "SELECT ".$os_taxTable.", description, go_ids FROM ".$view." WHERE seqId='".$seqId."'";
												#print "Executing the query: $query\n";
												$select_handle = $dbh->prepare_cached ($query);
												die "ERROR [$?]: can't prepare queries; aborting: ?\n" unless defined $select_handle;
													$select_handle->execute() or die $select_handle->errstr;
												@res = $select_handle->fetchrow_array;
												
												#Here print a warning if there is no results because there should always be organism/taxonomy and description 
												#for a protein
												if (scalar(@res) == 0){
													print LOG $seqId."\n";
													$discordantData++;
												 }#else{print "Query: $qName. Resulting GO ids from the SELECT statement: ".$res[2]."\n";}
										 }
									}
									
									#NEW CODE TO GET GO TERMS ASSOCIATED TO DOMAINS FROM THE VIEW
									#Extracts information only if $database is defined. (rpsblast)	
									if ((defined $database) and ($algorithmType eq "rpstblastn") ){
										#fetch informations from the database
										my $query = "SELECT go_ids FROM ".$view." WHERE cdName='".$seqId."'";
										#print "Executing SQL: $query";
										$select_handle = $dbh->prepare_cached ($query);
										die "ERROR [$?]: can't prepare queries; aborting: ?\n" unless defined $select_handle;
										$select_handle->execute() or die $select_handle->errstr;
										@res = $select_handle->fetchrow_array;
										#if (scalar(@res) > 0){
											#print "Result for $query: ".$res[0]." \n";
										 #}
									}
									
									#############################
									#Iterate on hsp with Bioperl
									while( my $hsp = $hit->next_hsp ) {
																		 
										#calculate query and hit coverage
										#HSP	 length('total')	  length of HSP (including gaps)
										#HSP	 length('hit')	  length of hit participating in alignment minus gaps
										#HSP	 length('query')		 length of query participating in alignment minus gaps
										#Calculate coverage. It will be displayed..  
										#The Query coverage is calculate as the length of the HSP 
										#that aligns on the query divided by the length of the query
										#multiplied per 100 to obtain a percentage  									                                      
										my $hspQueryCoverage = ($hsp->length('query')/$result->query_length)*100;
										#The Hit coverage is calculated as the length of the HSP 
										#that aligns on the HIT divided by the length of the HIT
										#multiplied per 100 to obtain a percentage   									 
										my $hspHitCoverage = ($hsp->length('hit')/$hit->length)*100;
										
										#Here we build the hash to have more fast informations extraction (only for BLASTX and BLASTN)		
										#ATTENTION (FUTURE VERSION): I could try to call some routine of BIOPERL to ask for the best scored hsp to insert... 	
										#and this control is also later in the GO and Descrs	 
										if( ($hspId == 1) and ($hitId == 1)){#Only the first result of the first hit has to be taken
											if ( ($algorithmType eq "blastx")){
											 # print "Algorithm is blastx!\n";
												$outHash->{$qName}->{$uniprotDBType}->{"name"} = $seqId;
												$outHash->{$qName}->{$uniprotDBType}->{"hspl"} = $hsp->length('query');#Length of the HSP without GAPS
												$outHash->{$qName}->{$uniprotDBType}->{"hitl"} = $hit->length;#Length of the HIT
												$outHash->{$qName}->{$uniprotDBType}->{'eval'} = $hsp->evalue;
												$outHash->{$qName}->{$uniprotDBType}->{"strand"} = $hsp->strand('query');#Strand of the query match
												$outHash->{$qName}->{$uniprotDBType}->{"frame"} = $hsp->query->frame;
												$outHash->{$qName}->{$uniprotDBType}->{"qc"} = $hspQueryCoverage;
												$outHash->{$qName}->{$uniprotDBType}->{"hc"} = $hspHitCoverage;
												 #print "name: ".$outHash->{$qName}->{$uniprotDBType}->{"name"}."\n";
											}				 
										}
										if ($algorithmType eq "blastn"){
											$outHash->{$qName}->{$algorithmType}->{"name"} = $seqId;
											$outHash->{$qName}->{$algorithmType}->{'eval'} = $hsp->evalue;
											$outHash->{$qName}->{$algorithmType}->{"hspl"} = $hsp->length('query');#Length of the HSP without GAPS
											$outHash->{$qName}->{$algorithmType}->{"strand"} = $hsp->strand('query');
											$outHash->{$qName}->{$algorithmType}->{"frame"} = $hsp->query->frame;
										}
													
										#Here we check the goodness of the HSP
										if( $hsp->length('total') > 1 ) {
											if ( $hsp->percent_identity >= 1 ) {
												
												#GFF PRINTING
												if ($configHash->{'printGFFOutput'} eq 'YES') {
													#match-part is a SOFA term used to describe a part of a match, for example an hsp from blast
													print GFF_FILE $qName."\t".$algorithmType."\t"."match_part\t";
													#my ($query_beg, $query_end) = $hsp->range('query');
													#my ($hit_beg, $hit_end) = $hsp->range('hit');
																	
													print GFF_FILE $hsp->start('query')."\t".#HSP START respect to the Query
																					$hsp->end('query')."\t";#END HSP respect to the Query
													printf GFF_FILE "%.2e\t",$hsp->evalue;#HSP SCORE
													
													#QUERY STRAND
													if ($hsp->strand('query') == "1"){
														print GFF_FILE "+\t";
													}else{print GFF_FILE "-\t";}
													
													if(defined $hsp->query->frame){
														print GFF_FILE $hsp->query->frame."\t";#QUERY FRAME
													}else{print GFF_FILE "\t";	}

													print GFF_FILE "Target=".$seqId." ".
													$hsp->start('query')." ".#HSP START respect to the Query
													$hsp->end('query').";".#HSP END respect to the Query
													"hin=".$hitId.";".#HIT NUMBER
													"hsn=".$hspId.";"; #HSP NUMBER
													if ($uniprotDBType ne 'none'){	print GFF_FILE "db=".$uniprotDBType.";";}
																					
													print GFF_FILE "hs=".$hsp->query->frame.";". #HSP STRAND
													"hf=".$hsp->hit->frame.";". #HIT FRAME
													"hsl=".$hsp->length('query').";". #HSP LENGTH respect to the Query
													"hil=".$hit->length.";". #HIT LENGTH
													"ql=".$result->query_length.";"; #QUERY LENGTH
													printf GFF_FILE "pi=%.2f;",$hsp->percent_identity; #PERCENTAGE IDENTITY HSP
													printf GFF_FILE "qc=%.2f;",$hspQueryCoverage; #QUERY COVERAGE of HSP
													printf GFF_FILE "hc=%.2f",$hspHitCoverage; #HIT COVERAGE of HSP
												}

																									
												#Extracts information from the UNIPROT database but only if $database is defined. (only BLASTX)
												if (defined $database and (scalar(@res) > 0)){
													if(defined($res[1]) and ($algorithmType eq "blastx")){
														my $description = $res[1];
														#We clean the description cancelling every strange character
														$description =~ s/[^A-Za-z0-9\-\(\)\s\_\/\.\']/ /g;
														if ($configHash->{'printGFFOutput'} eq 'YES') {
															print GFF_FILE ";d=".$description;
														}
														#only for the first result from the first HSP is taken
														if(($hspId == 1) and ($hitId == 1) and ($algorithmType eq "blastx")){
															$outHash->{$qName}->{$uniprotDBType}->{"d"} = $description;
														} 
													}
													
													if(defined($res[0]) and ($algorithmType eq "blastx") ){
														my $organism =$res[0];
														$organism =~ s/[^A-Za-z0-9\-\(\)\s\_\/\.\']/ /g;
														if ($configHash->{'printGFFOutput'} eq 'YES') {
															print GFF_FILE ";".$os_tax."=".$organism;
														}
														#only for the first result from the first HSP is taken
														if(($hspId == 1) and ($hitId == 1) ){
															$outHash->{$qName}->{$uniprotDBType}->{$os_tax} = $organism;
														} 
													}elsif ($algorithmType eq "blastx"){
														print "WARNING: organism for $qName has not been found in the database. ".
															"There's something wrong in the db writing or in the UNIPROT fasta file used as database. Check it!\n";
														}										
													
													if ( defined($res[2]) and ($algorithmType eq "blastx")){
														my $goIds = $res[2];

														#GFFPRINTING
														if ($configHash->{'printGFFOutput'} eq 'YES') {														
															#Conversion of numeric ids in GO terms using an SQL query 
															my $gffGOString = '';
															my @goTerms = split (/,/, $goIds);
															my @goRes = ();
															#FOR EACH GO TERM ASSIGNED TO THIS FEATURE
															#fetches informations from the uniprot database
															foreach my $goTerm (@goTerms){
																 my $go_select_handle = $dbh->prepare ("SELECT go_id FROM goTable WHERE id='".$goTerm."'");
																					 die "ERROR [$?]: can't prepare queries; aborting: ?"
																						unless defined $go_select_handle;
																 $go_select_handle->execute() or die $go_select_handle->errstr;
																 @goRes = $go_select_handle->fetchrow_array;	
																 if (scalar(@goRes) > 0){
																		$gffGOString.= $goRes[0].",";
																 }
															}	
															 chop($gffGOString);#remove the last ','
															 print GFF_FILE ";go=$gffGOString\n";
													 }
														
														#I save in the hash only the go terms for the first result
														if(($hspId == 1) and ($hitId == 1) and ($algorithmType eq "blastx") and
														($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both') ){
															$outHash->{$qName}->{$uniprotDBType}->{"go"} = $goIds;
															#print "GO Ids in the hash: ".$outHash->{$qName}->{$uniprotDBType}->{"go"}."\n";
														} 
													}#else{print GFF_FILE "\n";}
													
													#IF GO terms should be associated to domains and is defined the res[0]
													if ((defined $res[0]) and ($algorithmType eq "rpstblastn") ){
														my $goIds = $res[0];
														
														if ($configHash->{'printGFFOutput'} eq 'YES') {
															#GO terms associated to domains are always written in the gff
															#Getting the GO terms from the database to insert in the GFF file
															my $gffGOString = '';
															my @goTerms = split (/,/, $goIds);
															my @goRes = ();
															
															#FOR EACH GO TERM ASSIGNED TO THIS FEATURE
															#fetches informations from the uniprot database
															foreach my $goTerm (@goTerms){
																 my $go_select_handle = $dbh->prepare ("SELECT go_id FROM goTable WHERE id='".$goTerm."'");
																					 die "ERROR [$?]: can't prepare queries; aborting: ?"
																						unless defined $go_select_handle;
																 $go_select_handle->execute()
																					or die $go_select_handle->errstr;
																 @goRes = $go_select_handle->fetchrow_array;	
																 if (scalar(@goRes) > 0){
																		$gffGOString.= $goRes[0].",";
																 }
															}	
															chomp($gffGOString);
															print GFF_FILE ";go=".$gffGOString;
														}
														
														#We write in the hash all the GO terms associated to the HSPs separated with a symbol
														#only if the user chose to get them
														if (($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both') ){
															#$goIds =~ s/;/ /g;
															$outHash->{$qName}->{$algorithmType}->{"go"} .= $goIds.$separator;
														}
													
													} 
												}												
												if ( ($algorithmType eq "rpstblastn") or ($algorithmType eq "blastn") ) { #Otherwise it takes the hit->description (BLASTN)
													my $alternativeDescr = $hit->description;
													$alternativeDescr =~ s/[^A-Za-z0-9\-\(\)\s\_\/\.\']/ /g;					
													#GFF PRINTING
													if ($configHash->{'printGFFOutput'} eq 'YES') { 
														printf GFF_FILE ";d=".$alternativeDescr."\n";
													}
													if ( $algorithmType eq "blastn" ){ 
														$outHash->{$qName}->{$algorithmType}->{'d'} = $alternativeDescr;
													}
													#if ( $algorithmType eq "blastx" ){ 
														#$outHash->{$qName}->{$uniprotDBType}->{'d'} = $alternativeDescr;
													#}
												}
											}#ENDIF
										}#ENDIF
										$hspId++;
									}#END OF LOOP ON HSP
									if ( (defined $database) and (defined $select_handle ) ) {$select_handle->finish;}
								 $hitId++;
								}#END OF LOOP ON HIT
								
								 #When the loop finishes we don't need the last separator for the GO terms associated to domains
								 if (($algorithmType eq "rpstblastn") and
									 (($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both'))
									and (defined $outHash->{$qName}->{$algorithmType})){
										$outHash->{$qName}->{$algorithmType}->{"go"} =~ s/(\Q$separator\E)$//;
								 }
							}
						}#ENDWHILE ON THE RESULTS
					}#END TRY
					catch Bio::Root::Exception with {
						my $err = shift; 
						print "A Bioperl exception occurred:\n$err\n";
					};
				}#else {print "WARNING $blastOut has not been found!\n";}#DEBUGCODE
					#We check here if there is discordant data between the blast output and data in database
				if ( $discordantData > 0){
						print "Annocript has found $discordantData $uniprotDBType genes in the file $blastOut which are not present in the database."
								." Maybe you are using a BLAST output with version different from that of database"
								." Check this or contact us in the forum because something in the header of UniProt database might be changed.\n"
								." SQL query which produced the fail:\n SELECT $os_taxTable, description, go_ids FROM $view WHERE "
								." seqId='XXXXX') \n You can find all the gene accession ids in the log file at: $logFile.\n";
						$discordantData = 0;#Re-initialize
				}
			}
			if ($configHash->{'printGFFOutput'} eq 'YES') {
				print "\nYour blast output ".$fileName." has been converted in a gff3 file named ".$gff3FileName."\n";
				if ( $algorithmType eq "rpstblastn" ){
						print "GFF3 file for $algorithmType is always created to access data information faster.\n";
				}
			}
		}
	}else{print "WARNING: There are no output files in $outputsDir. GFF file will not be created...\n";}
	print "\n";
	

	#Close the log file
	close(LOG);
	#Delete the log file if it is empty
  if ( -s ($logFile) == 0){
    print "No errors in the parsing. $logFile will be erased.\n";
    delete_file($logFile);
  }
  
	if (defined($database)) {
		$dbh->commit();
		$dbh->disconnect(); 
	}
	if ($configHash->{'printGFFOutput'} eq 'YES') {
		close(GFF_FILE);
	}
	#Restore the variable
	if (defined $overridePrintGFF) {
		if ( $overridePrintGFF eq 'printGFF' )	{
			$configHash->{'printGFFOutput'} = 'NO';
		}
	}
return $return;
}



=head2 saveOutHash

 Title : saveOutHash
 Usage : saveOutHash( -file => 'gff3_formatted_file',
					);

 Function: saves the hash contained in the file 'outHashFile'

 Returns : a gff3 file

=cut
sub saveOutHash {
	my $outHashFile = shift; 		
  
	open(BXHASH,"> $outHashFile");
	print BXHASH Dumper($outHash); 
	close(BXHASH);
}


=head2 gff3_2_db

 Title : gff3_2_db
 Usage : gff3_2_db( -file => 'gff3_formatted_file',
					);

 Function: creates a database with the gff3 specifications. It gives any error message

 Returns : a gff3 file

=cut
sub gff3_2_db{

	my $gff3FileName = shift;
	my $dbname = shift;

	$gff3FileName = $gffFolder."/".$gff3FileName;
	print "Loading gff3 file inside a db named ".$dbname."...";
  
  #Creates a database
	my $gffBD = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
					 -dsn => $dsnGff, #"dbi:mysql:$dbname:$host:$port",
					 -user => $mySqlUser,
					 -pass=> $mySqlPass,
					 -create => 0);

  #Creates a loader
	my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store => $gffBD,
							 -verbose => 1,
							 -fast => 1,
							 -index_subfeatures => 1);

  
	if( $loader->load($gff3FileName)>0){
		print "...DONE!\n";
	 }else {print "WARNING: Zero sequences have been loaded!! (with BLASTN output can happen!)\n"}
}


=head2 fasta_2_gffdb

 Title : fasta_2_gffdb
 Usage : fasta_2_gffdb( -file => 'gff3_formatted_file',
					);

 Function: Inserts sequences from a fasta file inside a GFF3 database

 Returns : a gff3 file

=cut
sub fasta_2_gffdb{
	my $fa_file = shift;
	my $dbname = shift;
	my $create = 0;
	$create = 1 if shift eq 'Y';

  $fa_file = $sessionFolder."/".$fa_file;
  
  my $fa_name = extract_name($fa_file,"1");
   #Allowed characters for mysql name are [0-9,a-z,A-Z$_] 
  #Here we convert - and . to _
  
  print "Creating db ".$dbname."\n"; # with user:".$mySqlUser." and pwd: ".$mySqlPass."\n";
  # DATA SOURCE NAME
   
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { 'RaiseError' => 1 } );
  my $rc = $dbh->func("createdb", $dbname, $host, $mySqlUser, $mySqlPass, 'admin') 
    or die "Unable to create $dbname db. Try to solve the problem or contact the system administrator. Annocript  will stop. \n";
  
  
  #Creates a Store object associated with the GFF3 db given
  my $gffDB = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
						 -user => $mySqlUser,
						 -pass => $mySqlPass,
						 -dsn => $dsnGff,#"dbi:mysql:$dbname",
						 -create => $create);
						 
  #Creates a SeqIO object from the fasta file
  my $seqio = Bio::SeqIO->new(-format => 'fasta', -file => $fa_file);
	
  print "Number of seqs in the fasta file: ";
  my $command =  "fgrep -c '>' $fa_file";
  if ( try_exec_command ($command) < 1){
    die "ERROR: an error occurred while executing $command \n";
  }
  print "\n";
  
  my $totSeqs = 0;
	#With a WHILE loop we create a series of features from each of the elements of the fasta file
	my @features;
	while(my $seq = $seqio->next_seq) {
	 my $id = $seq->id;
   
   #This
	 $gffDB->insert_sequence($id,$seq->seq);

	 my $feature = Bio::SeqFeature::Generic->new(-seq_id => $id,
							-start => 1,
							-end => $seq->length,
							-primary => 'toplevel',
							-source => $fa_name,
							-strand => 1,
							-attributes => { ID => $id } );
	 push(@features,$feature);
   $totSeqs++;
	}
  
  print "Number of seqs stored: $totSeqs\n";
	#Inserts features inside the GFF DB
	$gffDB->store(@features) or die "\n ERROR [$?]: loading sequences in $dbname: ?\n";
	print "Your fasta file ".$fa_file." has been stored in ".$dbname." database\n";
}



=head2 printSequencesID

 Title : printSequencesID
 Usage : printSequencesID( -dbname => 'gff3 db name',
					);

 Function: print all the ids of sequences presents in the gff db. It's only debug function
 Returns : nothing

=cut
sub printSequencesID{
 my $dbname = shift;

 my $gffDB = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
						 -user => $mySqlUser,
						 -pass => $mySqlPass,
						 -dsn => $dsnGff,#"dbi:mysql:$dbname",
						 -create => 0);
 my @ids = $gffDB->seq_ids();
 print "Sequence ids: \n";
 foreach my $id (@ids){
	print $id."\n";
 }
}

=head2 rpstblastnInfos

 Title : rpstblastnInfos
 Usage : rpstblastnInfos( -db => 'name of the GFF3 db',
				-id 
			 );

 Function: Takes features from a GFF3 database and print them in a way that we choosed to represent.
		Algorithm choose based on the overlap between two sequences: if the missing part is little, it takes the 
		sequence with the best score. Changing a threshold you can decide how much restrictive has
		to be the choice.
    
    The concept is that if the sequence A is contained in some other seq B1,B2,..,BN, then the sequence will be excluded if 
    the missing part is less than a threshold and if the score is worsen than the other sequences overlapped.

		PSEUDOCODE:

		FOREACH feature 
			exclude = NO
			foverlapped= extract_overlapped
			IF |foverlapped| = 0 THEN take(feature)
			ELSE
				FOREACH fov in foverlapped
					IF overlap_lenght(feature, fov) >= threshold AND 
						score(feature) > score(fov) THEN
						exclude = YES
			IF exclude = NO	THEN
				take (Feature)		

 Returns : an array @finalLine with names, scores, descs, start&End for each of the sequences taken

=cut
sub rpstblastnInfos{
	#INPUT VARIABLES
	my $gffDB = shift;
	my $id = shift;
  my $NCflag = shift;#flag to use to find out if the sequence is non-coding
  
  
	my @finalLine = ();#Initialization
	my $finalNames = "";
	my $finalScores = "";
	my $finalDescs = "";		
	my $finalStartEnd = "";


	my $maxOverlap = $configHash->{'rpstMaxOverlap'};
  
  my @excludedList = ();#Array of sequences excluded
  
	#Extract all the features corresponding to the seqid in input and rpstblastn
	my @ssFeatures = $gffDB->features(-seqid => $id,-type => ['match_part:rpstblastn']);
	if(@ssFeatures>0){
		#For each of this features control if they are contained in the others
    #These objects are of type: Bio::SeqFeature
		foreach my $ssFeature (@ssFeatures){
      #print "\n".$id."- Sequence: ".$ssFeature->display_name." st:".$ssFeature->start." end: ".$ssFeature->end."\n";#DEBUG CODE
			#print "overlapped with: \n";
			#We put in some variables the needed fields of the sequence actually in exam	 		
			my $start = $ssFeature->start;
			my $end = $ssFeature->end;
			my $name = $ssFeature->display_name;
			my $score = $ssFeature->score;
      my $ssComplName =  $name.'_'.$start.'_'.$end;
      
      #Computation starts only if the sequence has not previously excluded because some other sequence can represent it
			if ( !( grep {/$ssComplName/} @excludedList) ){
        my $exclude = 0;
      
        #Extracts the sequences overlapped with ssFeature
        my @overlapped = $gffDB->features(-seqid => $id,
                  -type => ['match_part:rpstblastn'], 
                  -start => $start,
                  -end => $end,
                  -range_type => 'overlaps');

        #If the sequence A is contained in some other seq B1,B2,..,BN
        if (scalar(@overlapped) > 1){ #A sequence is everytime contained inside itself
          my $bestScore = $ssFeature->score;
          #For each of the overlaps the sequence will be excluded depending from the threshold we 	
          #used and if the score is worsen than the other sequences overlapped				
          foreach my $c (@overlapped){
            my $cName = $c->display_name;
            my $cComplName = $c->display_name.'_'.$c->start.'_'.$c->end;
            #to avoid the case that the sequence has been previously excluded
            if ( !( grep {/$cComplName/} @excludedList)){
               #to avoid the case that the sequence is the same domain
              if (not ( ($cName eq $name) and ($start == $c->start) and ($end == $c->end) ) ){
                
                 # print "Sequence: ".$c->display_name." st:".$c->start." end: ".$c->end."\n";#DEBUGCODE
                  #print "overlap dim: ".getMissingPart($start,$end,$c->start,$c->end)."\n";#DEBUGCODE
                  #We need that at least once the sequence is excluded. In that case it has a good representative among the remaining seqs
                  if (missOrKeepRPS_interval($ssFeature, $c, $maxOverlap) == 1){
                     # print "Exluding $ssComplName using  $cComplName...\n";#DEBUGCODE
                     $exclude = 1;
                   }
              }
            }#else{ print "$cName (".$c->start.",".$c->end.") is in the excluded list - NOT TO EVALUATE\n"}#DEBUGCODE
          }#ENDFOREACH
        }#ENDIF
        if ( $exclude == 1){
                  push(@excludedList,$ssComplName);#The names are not unique so this is a way to store unique names
                  #print "$ssComplName has been excluded!\n";#DEBUGCODE
        }
        if ( $exclude == 0){
          $finalStartEnd = $finalStartEnd.$ssFeature->start."-".$ssFeature->end.$separator;
          $finalNames = $finalNames.$ssFeature->display_name.$separator; 
          $finalScores = $finalScores.$ssFeature->score.$separator;
          $finalDescs = $finalDescs.join(" ",$ssFeature->get_tag_values('d')).$separator;	
        }
      }#else {print "$name ($start : $end) is in the excluded list - NOT TO EVAL\n";}#DEBUGCODE
		}#ENDFOREACH
    
	#Erase now the separator at the end of each string
	$finalStartEnd =~ s/(\Q$separator\E)$//;
	$finalNames =~ s/(\Q$separator\E)$//;
	$finalScores =~ s/(\Q$separator\E)$//;
	$finalDescs =~ s/(\Q$separator\E)$//;
	push (@finalLine,$finalNames,$finalStartEnd,$finalScores,$finalDescs);	
	}#ENDIF
	else{	
    push(@finalLine,"-","-","-","-");	
    $$NCflag++;#if there are no results the sequence could be non-coding
  }
	
	#print "cd desc: ".$cdDesc."\n";
	return @finalLine;
}


=head2 rpstblastnAll

 Title : rpstblastnAll
 Usage : rpstblastnAll( -db => 'name of the GFF3 db',
				-id 
			 );

 Function: Takes features from a GFF3 database and print them in a way that we choosed to represent.
            Every result is taken but the description will be written only once if the id is repeated.
            
 Returns : an array @finalLine with names, scores, descs, start&End for each of the sequences taken

=cut
sub rpstblastnAll{
	#INPUT VARIABLES
	my $gffDB = shift;
	my $id = shift;
  my $NCflag = shift;#flag to use to find out if the sequence is non-coding
  
  
	my @finalLine = ();#Initialization
	my $finalNames = "";
	my $finalScores = "";
	my $finalDescs = "";		
	my $finalStartEnd = "";
  
  my @descsIn = ();#Array of descriptions inserted
  
	#Extract all the features corresponding to the seqid in input and rpstblastn
	my @ssFeatures = $gffDB->features(-seqid => $id,-type => ['match_part:rpstblastn']);
	if(@ssFeatures>0){
		#For each of this features control if they are contained in the others
		foreach my $ssFeature (@ssFeatures){

      #print "\n".$id."- Sequence: ".$ssFeature->display_name." st:".$ssFeature->start." end: ".$ssFeature->end."\n";#DEBUG CODE
			#print "overlapped with: \n";
			#We put in some variables the needed fields of the sequence actually in exam	 		
			my $start = $ssFeature->start;
			my $end = $ssFeature->end;
			my $name = $ssFeature->display_name;
			my $score = $ssFeature->score;
      #my $ssComplName =  $name.'_'.$start.'_'.$end;
      
      my @description = $ssFeature->get_tag_values('d');
      $finalStartEnd = $finalStartEnd.$ssFeature->start."-".$ssFeature->end.$separator;
      $finalNames = $finalNames.(split(" ",$description[0]))[0].$separator; 
      $finalScores = $finalScores.$ssFeature->score.$separator;
        
      #If the sequence description has not been already inserted, do it
			if ( !( grep {/$name/} @descsIn) ){
        $finalDescs = $finalDescs.join(" ",@description).$separator;	
        push(@descsIn,$name);
      }
  
		}#ENDFOREACH
    
	#Erase now the separator at the end of each string
	$finalStartEnd =~ s/(\Q$separator\E)$//;
	$finalNames =~ s/(\Q$separator\E)$//;
	$finalScores =~ s/(\Q$separator\E)$//;
	$finalDescs =~ s/(\Q$separator\E)$//;
	push (@finalLine,$finalNames,$finalStartEnd,$finalScores,$finalDescs);	
	}#ENDIF
	else{	
    push(@finalLine,"-","-","-","-");	
    $$NCflag++;#if there are no results the sequence could be non-coding
  }
	
	#print "cd desc: ".$cdDesc."\n";
	return @finalLine;
}


=head2 blastnInfo

 Title : blastnInfo
 Usage : blastnInfo( -outHash => 'reference to the hash with data',
			-id => 'id of the sequence',
			 );

 Function: Takes features from an hash and print them in an array in order to print on the final table.

 Returns : an array with the data read for only the sequence given by 'id'

=cut
sub blastnInfo{
	#INPUT VARIABLES
	my $outHash = shift;
	my $id = shift;
	my $NCflag = shift;
	my @finalLine = ();#Initialization
	
	#print "Extracting blastn informations...\n";
	if(defined $outHash->{$id}->{'blastn'}){
		if(defined $outHash->{$id}->{'blastn'}->{'name'}) {
			push(@finalLine,$outHash->{$id}->{'blastn'}->{'name'});
    }else{push(@finalLine,"-");}
			
    if(defined $outHash->{$id}->{'blastn'}->{'eval'}) {
			push(@finalLine,$outHash->{$id}->{'blastn'}->{'eval'});
    }else{push(@finalLine,"-");}
			
		if(defined $outHash->{$id}->{'blastn'}->{'d'}) {
			push(@finalLine,$outHash->{$id}->{'blastn'}->{'d'});
    }else{push(@finalLine,"-");}
	}else {
    push(@finalLine,"-","-","-");
     $$NCflag++;#if there are no results the sequence could be non-coding
  }
	#print "cd desc: ".$cdDesc."\n";
	return @finalLine;
}


=head2 blastxInfosUniprot

 Title : blastxInfosUniprot
 Usage : blastxInfosUniprot( -outHash => 'reference to the hash',
			-id => 'identifier of the query in input',
      -uniprotDBType => 'the uniprot db type: sp or tr',
      -putOS => 'a flag: if it is OS_YES the organism will be print',
      -protDB => 'reference to database proteins'
			 );


 Function: Takes features from an hash and print them in an array in the order to print on the final table.

 Returns : an array with the data for only the sequence given by 'id'

=cut
sub blastxInfosUniprot{

  #INPUT VARIABLES
  my $outHash = shift;
  my $id = shift;
  my $uniprotDBType =shift;
  my $protDB = shift;
  
  #my $os = undef;
  my $uniprotViewEnzyme = 'uniprotkbViewEnzyme';
  my $uniprotViewPathways = 'uniprotkbViewPathways';
  
  #Evaluate time
  my $partTime = time;
  my $partDuration = undef;
	
	#print "Entering blastxinfotrsp\n";
  my @finalLine = ();#Initialization
  if(defined $outHash->{$id}->{$uniprotDBType}){
		#print "Printing information for $uniprotDBType\n";#DEBUGCODE
    #hsp Name
    if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
      push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'name'}); 
    }else{push(@finalLine, "-");}
				 	
		#prints hsp lenght if exists
    if(defined $outHash->{$id}->{$uniprotDBType}->{'hspl'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hspl'});
    }else{push(@finalLine, "-");}
				 
    #hsp evalue
    if(defined $outHash->{$id}->{$uniprotDBType}->{'eval'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'eval'});
    }else{push(@finalLine, "-");}
    
    #prints hit lenght if exists
    if(defined $outHash->{$id}->{$uniprotDBType}->{'hitl'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hitl'});
    }else{push(@finalLine, "-");}
    
    #print query coverage if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'qc'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'qc'});
		}else{push(@finalLine, "-");}
    
     #print hit coverage if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'hc'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'hc'});
		}else{push(@finalLine, "-");}
    
     #print strand if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'strand'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'strand'});
		}else{push(@finalLine, "-");}
		
		#print frame if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'frame'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'frame'});
		}else{push(@finalLine, "-");}
    		 
		#print description if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'d'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'d'});
		}else{push(@finalLine, "-");}
    
	}else{
		#If there is no output for this software it writes a set of '-' to indicate the absence
		push(@finalLine,"-","-","-","-","-","-","-","-","-");
	} 
	
  #The organism is written
	#Printing organism 
	if(defined $outHash->{$id}->{$uniprotDBType}->{$os_tax}){
		push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{$os_tax});
	}else{push(@finalLine, "-");}
		
	#IF Db is SwissProt add informations about the HSP using feature->display_name to query the protein db
	if( $uniprotDBType eq "sp" ){
	 if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
      #print "ID: $id dbtype $uniprotDBType - The name is: ".$outHash->{$id}->{$uniprotDBType}->{'name'}."\n";#DEBUGCODE
			push(@finalLine, enzymeInfo($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotViewEnzyme));
			#$partDuration = time - $partTime;
			#print "Enzyme occurs ".$partDuration." seconds\n";
			#$partTime = time;
			#FOR SINGLE COLUMN WITH PATHWAYS - COMMENTED
			#push(@finalLine, pathwaysInfo($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotViewPathways));
			#FOR MULTIPLE COLUMNS
			push(@finalLine, pathwaysInfoToCols($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotViewPathways));
			#$partDuration = time - $partTime;
			#print "Pathways occurs ".$partDuration." seconds\n";
			#$partTime = time;
			
	 }else {push(@finalLine,"-","-","-","-","-");}
 }	 
	return (\@finalLine);
}
 
=head2 blastxInfosUnirefSP

 Title : blastxInfosUnirefSP
 Usage : blastxInfosUnirefSP( -outHash => 'reference to the hash',
                    -id => 'identifier of the query in input',
                    -protDB => 'reference to database proteins'
                    );


 Function: Takes features from an hash and print them in an array in the order to print on the final table.
            This function has to be called when the Uniref database is used and fetches both the results 
            for SwissProt and Uniref.
            
 Returns : an array with the data for only the sequence given by 'id'

=cut
sub blastxInfosUnirefSP{

  #INPUT VARIABLES
  my $outHash = shift;#Hash where all the informations are contained
  my $id = shift;#Id of the sequence
  my $protDB = shift;#Reference to the db
  my $NCflag = shift;#flag to use to find out if the sequence is non-coding
  my $uniprotDBType = shift;#The db type

  my $os = undef; 
  my $uniprotViewEnzyme;
  my $uniprotWiewPathways;
  
  #Evaluate time
  my $partTime = time;
  my $partDuration = undef;

  #Let's change the view to use for access the enzyme corresponding
  #if ( $uniprotDBType eq 'uniref'){
	#$uniprotViewEnzyme = 'unirefViewEnzyme';
  #} else{
  $uniprotViewEnzyme = 'uniprotkbViewEnzyme';
  $uniprotWiewPathways = 'uniprotkbViewPathways';
  #}

	#print "Entering blastxinfotrsp\n";
  my @finalLine = ();#Initialization
  if(defined $outHash->{$id}->{$uniprotDBType}){
    #hsp Name
    if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
      push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'name'}); 
    }else{push(@finalLine, "-");}
				 	
    #prints hsp lenght if exists
    if(defined $outHash->{$id}->{$uniprotDBType}->{'hspl'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hspl'});
    }else{push(@finalLine, "-");}
				 
    #hsp evalue
    if(defined $outHash->{$id}->{$uniprotDBType}->{'eval'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'eval'});
    }else{push(@finalLine, "-");}
    
    #prints hit lenght if exists
    if(defined $outHash->{$id}->{$uniprotDBType}->{'hitl'}){
      push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hitl'});
    }else{push(@finalLine, "-");}
    
    #print query coverage if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'qc'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'qc'});
		}else{push(@finalLine, "-");}
    
     #print hit coverage if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'hc'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'hc'});
		}else{push(@finalLine, "-");}
    
    #print strand if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'strand'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'strand'});
		}else{push(@finalLine, "-");}
   
   	#print frame if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'frame'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'frame'});
		}else{push(@finalLine, "-");}
		 		 
		#print description if exist
		if(defined $outHash->{$id}->{$uniprotDBType}->{'d'}){
			push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'d'});
		}else{push(@finalLine, "-");}
    
	}else{
    #If there is no output for this software it writes a set of '-' to indicate the absence
    push(@finalLine,"-","-","-","-","-","-","-","-","-");
    $$NCflag++;#if there are no results the sequence could be non-coding
	} 
  
  #Printing taxonomy only for Uniref
  #if ( $uniprotDBType eq 'uniref'){
    if(defined $outHash->{$id}->{$uniprotDBType}->{$os_tax}){
       $os= $outHash->{$id}->{$uniprotDBType}->{$os_tax};
    }else{$os = '-'}
    push(@finalLine, $os);
  #}
  
  #Enzyme and Pathways information only for SwissProt
  if ( $uniprotDBType eq 'sp'){
    if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
      push(@finalLine, enzymeInfo($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotViewEnzyme));
      #$partDuration = time - $partTime;
      #print "Enzyme occurs ".$partDuration." seconds\n";
      #$partTime = time;
      
      #FOR SINGLE COLUMN WITH PATHWAYS - COMMENTED
      #push(@finalLine, pathwaysInfoToCols($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotWiewPathways));
      #THREE COLUMNS FOR PATHWAYS
      push(@finalLine, pathwaysInfoToCols($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$uniprotWiewPathways));
      
      #$partDuration = time - $partTime;
      #print "Pathways occurs ".$partDuration." seconds\n";
      #$partTime = time;
    }else {push(@finalLine,"-","-","-","-","-");}
  }
	return @finalLine;
}


#=head2 blastxInfosUniref

 #Title : blastxInfosUniref
 #Usage : blastxInfosUniref( -outHash => 'reference to the hash',
                    #-id => 'identifier of the query in input',
                    #-protDB => 'reference to database proteins'
                    #);


 #Function: Takes features from an hash and print them in an array in the order to print on the final table.

 #Returns : an array with the data for only the sequence given by 'id'

#=cut
#sub blastxInfosUniref{

  ##INPUT VARIABLES
  #my $outHash = shift;#Hash where all the informations are contained
  #my $id = shift;#Id of the sequence
  #my $protDB = shift;#Reference to the db
  #my $NCflag = shift;#flag to use to find out if the sequence is non-coding
  
  #my $uniprotDBType = 'uniref';
  #my $unirefViewEnzyme = 'unirefViewEnzyme';
  #my $os = undef; 
  
  ##Evaluate time
  #my $partTime = time;
  #my $partDuration = undef;

	##print "Entering blastxinfotrsp\n";
  #my @finalLine = ();#Initialization
  #if(defined $outHash->{$id}->{$uniprotDBType}){
    ##hsp Name
    #if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
      #push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'name'}); 
    #}else{push(@finalLine, "-");}
				 	
    ##prints hsp lenght if exists
    #if(defined $outHash->{$id}->{$uniprotDBType}->{'hspl'}){
      #push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hspl'});
    #}else{push(@finalLine, "-");}
				 
    ##hsp evalue
    #if(defined $outHash->{$id}->{$uniprotDBType}->{'eval'}){
      #push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'eval'});
    #}else{push(@finalLine, "-");}
    
    ##prints hit lenght if exists
    #if(defined $outHash->{$id}->{$uniprotDBType}->{'hitl'}){
      #push(@finalLine, $outHash->{$id}->{$uniprotDBType}->{'hitl'});
    #}else{push(@finalLine, "-");}
    
    ##print query coverage if exist
		#if(defined $outHash->{$id}->{$uniprotDBType}->{'qc'}){
			#push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'qc'});
		#}else{push(@finalLine, "-");}
    
     ##print hit coverage if exist
		#if(defined $outHash->{$id}->{$uniprotDBType}->{'hc'}){
			#push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'hc'});
		#}else{push(@finalLine, "-");}
    
     ##print strand if exist
		#if(defined $outHash->{$id}->{$uniprotDBType}->{'strand'}){
			#push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'strand'});
		#}else{push(@finalLine, "-");}
    		 
		##print description if exist
		#if(defined $outHash->{$id}->{$uniprotDBType}->{'d'}){
			#push(@finalLine,$outHash->{$id}->{$uniprotDBType}->{'d'});
		#}else{push(@finalLine, "-");}
    
	#}else{
    ##If there is no output for this software it writes a set of '-' to indicate the absence
    #push(@finalLine,"-","-","-","-","-","-","-","-");
    #$$NCflag++;#if there are no results the sequence could be a long non-coding
	#} 

	##Printing taxonomy 
	#if(defined $outHash->{$id}->{$uniprotDBType}->{$os_tax}){
		 #$os= $outHash->{$id}->{$uniprotDBType}->{$os_tax};
	#}else{$os = '-'}
  #push(@finalLine, $os);
	##Enzyme information
	 #if(defined $outHash->{$id}->{$uniprotDBType}->{'name'}){
			#push(@finalLine, enzymeInfo($outHash->{$id}->{$uniprotDBType}->{'name'},$protDB,$unirefViewEnzyme));
			##$partDuration = time - $partTime;
			##print "Enzyme occurs ".$partDuration." seconds\n";
			##$partTime = time;
			
	 #}else {push(@finalLine,"-","-");}
		 	 
	#return @finalLine;
#}


=head2 blastxGOInfos

 Title : blastxGOInfos
 Usage : blastxGOInfos( -outHash => 'an hash where all the info are written',
                        -dbhProt => 'reference to the db',
                        -id => 'identifier of the query in input'
                        - uniprotDBType => 'SP, swiss-prot; TR, trembl; UNIREF, uniref'
                       );

 Function: Takes GO terms from the main hash and print them in an array in the order to print on the final table.
            This task is computed in two different ways according to the uniprotDBType: uniprot_kb or uniref.

 Returns : an array with the data for only the sequence given by 'id'

=cut
sub blastxGOInfos{
	my $outHash = shift;
	my $id =shift;
	my $uniprotDBType =shift;
	my $dbhProt = shift;
	 
	#WHILE 
	my @finalLine = ();#Initialization
	
	if(defined $outHash->{$id}->{$uniprotDBType}){
		 #Extract values from the TAG "go" - BIO::SeqFeature::Generic
		 if(defined $outHash->{$id}->{$uniprotDBType}->{'go'}){
	 		#print "GO Terms: ".$values[0]."\n";
			my $bp = ""; my $mf = ""; my $cc = ""; 
			my $bpDef = ""; my $mfDef = ""; my $ccDef = "";
      
		 	my @goTerms = split (/,/, $outHash->{$id}->{$uniprotDBType}->{'go'});
			my @res = ();
      
			#FOR EACH GO TERM ASSIGNED TO THIS FEATURE
			#fetches informations from the uniprot database
			foreach my $goTerm (@goTerms){
			 my $select_handle = $dbhProt->prepare ("SELECT division, definition, go_id FROM goTable WHERE id='".$goTerm."'");
								 die "ERROR [$?]: can't prepare queries; aborting: ?"
									unless defined $select_handle;
			 $select_handle->execute()
								or die $select_handle->errstr;
			 @res = $select_handle->fetchrow_array;
  
			 if (scalar(@res) > 0){
					 # print "Division: ".$res[0]." Def: ".$res[1]."\n";#DEBUG CODE
					 #if arguments are not null insert them in variables
						 if ($res[0] eq "P") { $bp=$bp.$res[2].$separator; $bpDef=$bpDef.$res[1].$separator}
						 if ($res[0] eq "F") { $mf=$mf.$res[2].$separator; $mfDef=$mfDef.$res[1].$separator}
						 if ($res[0] eq "C") { $cc=$cc.$res[2].$separator; $ccDef=$ccDef.$res[1].$separator}	
						 #else	{ print "value not provided\n" }
				}#ENDIF 
						 					 	
			}#ENDFOREACH
			if (scalar(@res) > 0){
				 #Erasing all the separators at the end of the string
				 $bp =~ s/(\Q$separator\E)$//;
				 $mf =~ s/(\Q$separator\E)$//;
				 $cc =~ s/(\Q$separator\E)$//;
				 $bpDef =~ s/(\Q$separator\E)$//;
				 $mfDef =~ s/(\Q$separator\E)$//;
				 $ccDef =~ s/(\Q$separator\E)$//;
			}
				#print "\nbp: ".$bp." - mf: ".$mf." - cc: ".$cc."\n";
				#print "\nbpDef: ".$bpDef." - mfdef: ".$mfDef." - ccDef:".$ccDef."\n";
				#The previous variable are used to build strings for bp,mf,cc
		if ($bp ne "") {
          push(@finalLine, $bp);
        } else {push(@finalLine, "-");};
		if ($bpDef ne "") {
          push(@finalLine, $bpDef);
        } else {push(@finalLine, "-");};
		if ($mf ne "") {
          push(@finalLine, $mf);
        } else {push(@finalLine, "-");};
		if ($mfDef ne "") {
          push(@finalLine, $mfDef);
        } else {push(@finalLine, "-");};
		if ($cc ne "") {
          push(@finalLine, $cc);
        } else {push(@finalLine, "-");};
		if ($ccDef ne "") {
          push(@finalLine, $ccDef);
        } else {push(@finalLine, "-");};
						
		}else{ push(@finalLine, "-","-","-","-","-","-");}
  }else{ push(@finalLine, "-","-","-","-","-","-");}
	return @finalLine;
}
	
=head2 rpstblastnGOInfos

 Title : rpstblastnGOInfos
 Usage : rpstblastnGOInfos( -outHash => 'an hash where all the info are written',
                        -dbhProt => 'reference to the db',
                        -id => 'identifier of the query in input'
                        - uniprotDBType => 'SP, swiss-prot; TR, trembl; UNIREF, uniref'
                       );

 Function: Takes GO terms from the main hash and print them in an array in the order to print on the final table.
            This task is computed in two different ways according to the uniprotDBType: uniprot_kb or uniref.

 Returns : an array with the data for only the sequence given by 'id'

=cut
sub rpstblastnGOInfos{
	my $outHash = shift;
	my $id =shift;
	my $algorithmType =shift;
	my $dbhProt = shift;
	 
	#WHILE 
	my @finalLine = ();#Initialization
	my $semicolon = ';';
	
	#Here we have to separate the GO ids using commas and separator and then re-order them with the same approach
	if(defined $outHash->{$id}->{$algorithmType}){
		 #Extract values from the TAG "go" - BIO::SeqFeature::Generic
		 if(defined $outHash->{$id}->{$algorithmType}->{'go'}){
	 		#print "GO Terms: ".$values[0]."\n";
			my $bp = ""; my $mf = ""; my $cc = ""; 
			my $bpDef = ""; my $mfDef = ""; my $ccDef = "";
			
      #Array containing the GO terms inserted in the result
      my @goInserted = ();
			#We should make a loop that iterates on elements divided with the separator
			my @domainsGOs = split(/\Q$separator\E/, $outHash->{$id}->{$algorithmType}->{'go'});
		 	foreach my $domainGOs (@domainsGOs){
				#then another loop iterating on elements divided by comma. Those divided by comma
				my @goTerms = split (/,/, $domainGOs);
				my @res = ();
						  
				#FOR EACH GO TERM ASSIGNED TO THIS FEATURE
				#fetches informations from the uniprot database
				foreach my $goTerm (@goTerms){
					 my $select_handle = $dbhProt->prepare ("SELECT division, definition, go_id FROM goTable WHERE id='".$goTerm."'");
										 die "ERROR [$?]: can't prepare queries; aborting: ?"
											unless defined $select_handle;
					 $select_handle->execute() or die $select_handle->errstr;
					 @res = $select_handle->fetchrow_array;
					
					 
					 if (scalar(@res) > 0){
						 #If the go term as not been already inserted for this block of domains...
						 if ( scalar(grep {/\b$goTerm\b/} @goInserted) == 0){
							 # print "Division: ".$res[0]." Def: ".$res[1]."\n";#DEBUG CODE
							 #if arguments are not null insert them in variables
							 if ($res[0] eq "P") { $bp=$bp.$res[2].$separator; $bpDef=$bpDef.$res[1].$separator}
							 if ($res[0] eq "F") { $mf=$mf.$res[2].$separator; $mfDef=$mfDef.$res[1].$separator}
							 if ($res[0] eq "C") { $cc=$cc.$res[2].$separator; $ccDef=$ccDef.$res[1].$separator}		
							 push(@goInserted,$goTerm);#Insert this go term in the list of the already inserte				 
							}
					 }#ENDIF 
														
				}#ENDFOREACH
			}
			#Erasing all the separators at the end of the strings
			$bp =~ s/(\Q$separator\E)$//;
			$mf =~ s/(\Q$separator\E)$//;
			$cc =~ s/(\Q$separator\E)$//;
			$bpDef =~ s/(\Q$separator\E)$//;
		  $mfDef =~ s/(\Q$separator\E)$//;
			$ccDef =~ s/(\Q$separator\E)$//;
		
			#print "\nbp: ".$bp." - mf: ".$mf." - cc: ".$cc."\n";
			#print "\nbpDef: ".$bpDef." - mfdef: ".$mfDef." - ccDef:".$ccDef."\n";
			#The previous variable are used to build strings for bp,mf,cc
			if ($bp ne "") {
			  push(@finalLine, $bp);
			} else {push(@finalLine, "-");};
			if ($bpDef ne "") {
			  push(@finalLine, $bpDef);
			} else {push(@finalLine, "-");};
			if ($mf ne "") {
			  push(@finalLine, $mf);
			} else {push(@finalLine, "-");};
			if ($mfDef ne "") {
			  push(@finalLine, $mfDef);
			} else {push(@finalLine, "-");};
			if ($cc ne "") {
			  push(@finalLine, $cc);
			} else {push(@finalLine, "-");};
			if ($ccDef ne "") {
			  push(@finalLine, $ccDef);
			} else {push(@finalLine, "-");};
						
		}else{ push(@finalLine, "-","-","-","-","-","-");}
  }else{ push(@finalLine, "-","-","-","-","-","-");}
	return @finalLine;
}	

  
=head2 collectFastaNames

 Title : collectFastaNames
 Usage : collectFastaNames( -fastaFile => 'the fasta where to collect',
			 );

 Function: collects the identifiers from a FASTA file
 Returns : an array of identifiers

=cut
sub collectFastaNames{
  my $fastaFile = shift;
  
  my $fastaPath = $sessionFolder."/".$fastaFile;
  my @ids = ();
  open(FASTA,"<$fastaPath") or die "ERROR: Cannot open $fastaPath\n";
  
  while(my $line = <FASTA>){
      if( $line =~ /\>/){
          my $id = (split(" ", $line))[0];#Collects the name comprising '>'
          $id = substr($id,1);#Remove the '>' from the name
          push(@ids, $id);
      }
  }
  close(FASTA);
  return @ids;
}  
	

=head2 printAllInfo

 Title : printAllInfo
 Usage : printAllInfo( -dbname => 'name of the database',
			-fileTableName => 'name of the Annocript final output file',
			-outHashFile => 'the file where the hash is stored',
			 );

 Function: This is the main subroutine to write the complete output file. It will flow the various hashes prepared with the
          proper informations. It uses both an hash and a GFF database to extract them.
           It starts by loading the hash from file and accessing the protein database and GFF3 one. It iterate thorugh the 
            various sequences of the fasta file given and, for each of the, extract the corresponding results for all the analyses
          obtain.
          
          The EnzymeIDs and Pathways are given only for SwissProt.
          When the table is referred to Uniprot_KB database the GO term is associated with the lowest value between SP and TREMBL.
          When the table is referred to Uniref database the GO term is associated with only Uniref and not SP.
          Six conditions have to be satisfied from a sequence to be non coding. In this function we also insert this conditions and
          thanks to a variable we check if them all are satisfied.
          PAY ATTENTION:: REMEMBER TO CHANGE THE NUMBER OF CONDITIONS TO SATISFY IF YOU CHANGE THE CRITERIA!!!!
 
 Returns : a file tab separated with the datas

=cut
sub printAllInfo{
 my $dbname = shift;
 my $fileTableName = shift; 
 my $ridFileTable = shift;
 #my $storeBlastx = shift;
 my $outHashFile = shift;
 
 
 #Variables to see time occurrence
 my $partTime = time;
 my $partDuration = undef;
	
	
 my $gffDB = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
						 -user => $mySqlUser,
						 -pass => $mySqlPass,
						 -dsn => $dsnGff,#"dbi:mysql:$dbname;mysql_local_infile=1",
						 -create => 0);
 # PERL DBI CONNECT TO PROTEIN DB
 my $proteinDB = DBI->connect($dsn, $mySqlUser, $mySqlPass, { PrintError => 1, AutoCommit => 0 });
	
 
  #Loading the hash from a file if it hasn't been created
  #It does not contain the information for RPSBLAST!!!!
  #if (!(defined $outHash)){
		#$outHash = retrieve ($outHashFile);
  #} 
  
  print "Creating the final Annocript output file...\n";
  #FOR EACH SEQUENCE RECORDED IN THE GFF DB. I need this array to store the sequences
  #but I will not use it in the loop since I want to use the same loop if the user
  #has a filtered fasta file and wants an output only with that
  my @ids = $gffDB->seq_ids();
  
  #Instead, the loop we'll run using the names in the fasta file in input
  my @ids4PrintInfo = @tr_ids;#collectFastaNames($configHash->{'fastaSeqs'});
  print "(".scalar(@ids4PrintInfo)." sequences) Sequences count: ";
	
  $fileTableName = $outFolder."/".$fileTableName;
  $ridFileTable = $outFolder."/".$ridFileTable;
  
  #Opening a file where  results can be print
  open (OUT_FILE,">$fileTableName") || die ( "\nERROR [$?]: Cannot open file $fileTableName to write: ?" );
  open (OUT_RID_FILE,">$ridFileTable") || die ( "\nERROR [$?]: Cannot open file ridFileTable to write: ?" );
  #Printing on the file the header (ADD SOMETHING HERE IF YOU WANT TO ADD FIELDS!!

	
  #Number of the conditions to satisfy for the sequence to be non-coding
  #NB: IF YOU ARE ADDING SOME FIELD THAT INDICATES NONCODINGness THAN CHANGE THE VALUE ACCORDINGLY
  my $totConditions = 7;
   
  #my $counter = 1;	
  my $countThresh = 1000;
  my $id = '';
  my $sequence = '';
  my $seqLength = 0;
  my $chosenDB = "";#Variable used to select the database of proteins from the Hash (sp,tr,uniref)
  
  for (my $counter = -1; $counter < scalar(@ids4PrintInfo); $counter++){
		#foreach my $id (@ids4PrintInfo){
		#print "Counter of lines: $counter\n";
		#COMMON FEATURES ELEMENTS
		if ($counter >= 0){	
			$id = 	$ids4PrintInfo[$counter];
			#print "DEBUG CODE: Loop for write the output- Id = $id\n";
			$sequence = $gffDB->fetch_sequence( -seq_id => $id );
			$seqLength = length($sequence);
		}
		my $NCflag = 0;#This flag will be activated if some conditions (inserted in the subroutines) are verified

		#my @features = $gffDB->features( -seqid => $id );

		
		#Print the columns header
		if ($counter == -1){
			#print "Printing the header ".join("\t",@transcriptsCharHeader)."\n";
			print OUT_FILE join("\t",@transcriptsCharHeader);
			print OUT_FILE "\t";
			
			print OUT_RID_FILE join("\t",@transcriptsCharHeader);
			print OUT_RID_FILE "\t";	  
		}else{		
			print OUT_FILE $id."\t";
			print OUT_RID_FILE $id."\t";
			print OUT_FILE $seqLength."\t";
			print OUT_RID_FILE $seqLength."\t";
		}	
		
		#Depending from the number of sequences the current state will be printed
		my $totSequences = $configHash->{'totSequences'};
		if ( $totSequences < $countThresh){
			$countThresh = 100;
			if ( $totSequences < $countThresh){
			$countThresh = 10;
			if ( $totSequences < $countThresh){
				$countThresh = 1;
			}
			}
		}	
		#Prints the number of sequences evaluated
		if (($counter >= 0) and ($counter % $countThresh == 0) ){
			print $counter." - ";
		}
		#$counter++;#Counts the number of sequences evaluated
			
		#The first condition of the non codingness of a sequence is about its lenght
		if ( ($counter >= 0) and ($seqLength >= $configHash->{'NCSeqLength'}) ){
			$NCflag++;
		}
			
		
		#########################
		#Here we print informations related to SwissProt and UniRef or TrEMBL
		#Print the columns header
		if ($counter == -1){
			if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){
				$chosenDB = "tr";
				#print "Printing the header ".join("\t",@blastxSPTRHeader)."\n";
				print OUT_FILE join("\t",@blastxSPTRHeader);
				print OUT_FILE "\t";
				
				print OUT_RID_FILE join("\t",@blastxSPTRHeader);
				print OUT_RID_FILE "\t";	  
				
				#If the GO terms should be print for proteins, print the header now 
				if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
					#print "Printing the header ".join("\t",@protGOHeader)."\n";
					print OUT_FILE join("\t",@protGOHeader);
					print OUT_FILE "\t";
				
					print OUT_RID_FILE join("\t",@protGOHeader);
					print OUT_RID_FILE "\t";
			 }
			}  
			if ( $configHash->{'dbInUse'} eq 'uniref' ){
				 $chosenDB = "uniref";
					#Print the columns header
					print OUT_FILE join("\t",@blastxSPUfHeader);
					print OUT_FILE "\t";
					
					print OUT_RID_FILE join("\t",@blastxSPUfHeader);
					print OUT_RID_FILE "\t";	  
					
					#If the GO terms should be print for proteins, print the header now 
					if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
						#print "Printing the header ".join("\t",@protGOHeader)."\n";
						print OUT_FILE join("\t",@protGOHeader);
						print OUT_FILE "\t";
					
						print OUT_RID_FILE join("\t",@protGOHeader);
						print OUT_RID_FILE "\t";
				 }
			}
		}else{
			
			#Here I Print the best result (lowest evalue) among the HPSs found to be homologous 			
			#my $start = time();
			#Print a separated result for BLASTX on trembl and Swiss-prot
			#This includes the print of Enzyme infos
			my @blastxElements = printBlastProteinBestResults($outHash,$id,$proteinDB,\$NCflag,$chosenDB);
			print OUT_FILE join("\t",@blastxElements);
			print OUT_FILE "\t";

			print OUT_RID_FILE join("\t",@blastxElements);
			print OUT_RID_FILE "\t";
			#my $end = time();
			#printf(" blastx :%.2f\n", $end - $start);	 					
		}
			
		#Print the columns header
		if ($counter == -1){
			#print "Printing the header ".join("\t",@rpstblastnHeader)."\n";			
			print OUT_FILE join("\t",@rpstblastnHeader);
			print OUT_FILE "\t";
				
			print OUT_RID_FILE join("\t",@rpstblastnHeader);
			print OUT_RID_FILE "\t";	
			
			if ( ($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both')){
				print OUT_FILE join("\t",@domGOHeader);
				print OUT_FILE "\t";
			
				print OUT_RID_FILE join("\t",@domGOHeader);
				print OUT_RID_FILE "\t";
			}
		}else{
			#Print on a line informations retrieved from the rpsblastn output
			#print "Printing rpstblastn informations...";
			#print OUT_FILE join("\t",rpstblastnInfos($gffDB,$id,\$NCflag));#Function to separate intervals
			my @rpsElements =  rpstblastnAll($gffDB,$id,\$NCflag);
			print OUT_FILE join("\t",@rpsElements);
			print OUT_FILE "\t";
			
			print OUT_RID_FILE join("\t",@rpsElements);
			print OUT_RID_FILE "\t";
			#print "...DONE!";
		
			if ( ($configHash->{'goTermsAss'} eq 'domains') or ($configHash->{'goTermsAss'} eq 'both')){
				my @rpsGOElements = rpstblastnGOInfos($outHash,$id,'rpstblastn',$proteinDB);
				
				print OUT_FILE join("\t",@rpsGOElements);
				print OUT_FILE "\t";
				
				print OUT_RID_FILE join("\t",@rpsGOElements);
				print OUT_RID_FILE "\t";
				#print "...DONE!";
			}	
		}    

		#Print the columns header
		if ($counter == -1){
			print OUT_FILE join("\t",@blastnHeader);
			print OUT_FILE "\t";
				
			print OUT_RID_FILE join("\t",@blastnHeader);
			print OUT_RID_FILE "\t";	
		}else{
			my @blastnElements =  blastnInfo($outHash,$id,\$NCflag);
			#Print on a line informations retrieved from the blastn output
			#print "Printing blastn informations...";
			print OUT_FILE join("\t",@blastnElements);
			print OUT_FILE "\t";
			
			print OUT_RID_FILE join("\t",@blastnElements);
			print OUT_RID_FILE "\t";
			#print "...DONE!\n";
			#print "blastn-The NCflag is: $NCflag\n";		
		}		
		

		#Print the columns header
		if ($counter == -1){
			print OUT_FILE join("\t",@orfHeader);
			print OUT_FILE "\t";
				
			print OUT_RID_FILE join("\t",@orfHeader);
			print OUT_RID_FILE "\t";	
		}else{
			my @orfElements = orfInfo($outHash,$id,\$NCflag);
			#Print on a line informations retrieved from the dna2pep output
			#print "Printing ORF informations...";
			print OUT_FILE join("\t", @orfElements);
			print OUT_FILE "\t";
			
			print OUT_RID_FILE join("\t",@orfElements);
			print OUT_RID_FILE "\t";
			#print "...DONE!\n";	
		}

		#Print the columns header
		if ($counter == -1){
			print OUT_FILE join("\t",@ncHeader);
			print OUT_FILE "\t";
				
			print OUT_RID_FILE join("\t",@ncHeader);
			print OUT_RID_FILE "\t";	
		}else{
			my $ncElement;
			$ncElement = nonCodingInfo($outHash,$id,\$NCflag);
			#Print on a line informations retrieved from the lnc prediction output
			#print "Printing lnc prediction informations...";
			print OUT_FILE join("\t", $ncElement);
			print OUT_FILE "\t";
			
			print OUT_RID_FILE join("\t", $ncElement);
			print OUT_RID_FILE "\t";
			#print "...DONE!\n";  
		}

		#Print the columns header
		if ($counter == -1){
			print OUT_FILE join("\t",@heuristicHeader);
			print OUT_FILE "\t";
				
			print OUT_RID_FILE join("\t",@heuristicHeader);
			print OUT_RID_FILE "\n";	
		}else{
			#Prints our heuristical prediction of long non-coding (1:non-coding)
			if ($NCflag == $totConditions){ 
				#print "The sequence $id is NC $NCflag constrains are satisfied\n";
				print OUT_FILE '1';
				print OUT_RID_FILE '1';
			}else{
				print OUT_FILE '0';
				print OUT_RID_FILE '0';
				}
				print OUT_FILE "\t";
				print OUT_RID_FILE "\n";
		}
			
		#Print the columns header
		if ($counter == -1){
			print OUT_FILE join("\t",@seqHeader);
			print OUT_FILE "\n";
		}else{
			#Prints the sequence of all the query
			print OUT_FILE $sequence;
			print OUT_FILE "\n";
			}	
	}

	print "DONE! \nThe output has been written in ".$fileTableName.".\n";
	close(OUT_FILE);
  close(OUT_RID_FILE);
	$proteinDB->commit();
	$proteinDB->disconnect(); 
}


								
=head2 enzymeInfo

 Title : enzymeInfo
 Usage : enzymeInfo( -dbname => 'name of the GFF3 db',
			-use => 'username for dbname',
			-pass => 'password for dbname',
			-dbProt => 'the database uniprot created',
			-userProt => 'username for uniprot',
			-passProt => 'password for uniprot'
			 );

 Function: Returns an array with two elements that will be columns: one is the enzyme ids for a given protein and the other is the descriptions.
			Takes in input an uniprot id that is usually the HSP name.
 
 Returns : an array

=cut
sub enzymeInfo{

  #VARIABLES FOR USE OF THE UNIPROT DATABASE
	my $uniprotId = shift;
	my $dbh = shift;
	my $viewEnzyme = shift;
  
  my @res = ();
  #Accessing database to extract enzyme ids corresponding to a certain uniprot id
	my $select_handle = $dbh->prepare("SELECT enzyme_ids FROM $viewEnzyme WHERE seqId='".$uniprotId."'");
	print "WARNING! In enzymeInfo: Couldn't prepare first query; " unless defined $select_handle;
  $select_handle->execute()	or print $select_handle->errstr;
	@res = $select_handle->fetchrow_array;
							
	my @finalLine = ();#Initialization of the line that will be printed
	
	# print "uniprotid: ".$uniprotId."\n";
  if ( defined $res[0] ){
    #print "Enzyme_ids: ".$res[0]."\n"; #DEBUGCODE
    my $ids = "";
    my $descs = "";
    my @enzymeIds = split(",",$res[0]);#The enzyme ids can be a list separated by comma
    #but these are database ids, we must find the corresponding names in the table
    foreach my $enzymeId (@enzymeIds){
      #print "Doing the query to extract seqId...\n";
      #print "Second Enzyme query";
      my $select_handle2 = $dbh->prepare("SELECT enz_id, description FROM enzymeTable WHERE id='".$enzymeId."'");
      print "WARNING! In enzymeInfo: Couldn't prepare second query; "
      unless defined $select_handle2; 
      $select_handle2->execute() or print $select_handle2->errstr;
      my @res2 = $select_handle2->fetchrow_array;
      
      #from this second result we can create the final string to add
      if (scalar(@res2) > 0){
        #This IF is created because sometimes something has already been
        #written. If this is the case then separate with comma
        #otherwise it is the first element
				if ($ids eq ''){
          $ids = $res2[0];
        }else{$ids = $ids.$separator.$res2[0];}
				
        if ($descs eq ''){
          $descs = $res2[1];
        }else{$descs = $descs.$separator.$res2[1];}
      }
      #We need this else when the query gives an empty result. This means that there's no enzyme corresponding to the id
      else{
        print "ERROR in fetching informations from EnzymeTable in database. The id= $enzymeId did not correspond to anything.Check integrity of database!\n";
			 }
    }
    push ( @finalLine, $ids, $descs); 
		
  }else {push ( @finalLine,"-","-");}
	return @finalLine;
}
	
								
=head2 pathwaysInfo

 Title : pathwaysInfo
 Usage : pathwaysInfo( -dbname => 'name of the GFF3 db',
			-use => 'username for dbname',
			-pass => 'password for dbname',
			-dbProt => 'the database uniprot created',
			-userProt => 'username for uniprot',
			-passProt => 'password for uniprot'
			 );

 Function: Returns the descriptions of the pathways.
			Takes in input an uniprot id that is usually the HSP name.
 
 Returns : an array

=cut
sub pathwaysInfo{

  #VARIABLES FOR USE OF THE UNIPROT DATABASE
	my $uniprotId = shift;
	my $dbh = shift;
	my $viewPathway = shift;
  
  my @res = ();
  #Accessing database to extract enzyme ids corresponding to a certain uniprot id
	my $select_handle = $dbh->prepare("SELECT pw_ids FROM $viewPathway WHERE seqId='".$uniprotId."'");
	print "WARNING! In pathwaysInfo: Couldn't prepare first query; " unless defined $select_handle;
  $select_handle->execute()	or print $select_handle->errstr;
	@res = $select_handle->fetchrow_array;
							
	my @finalLine = ();#Initialization of the line that will be printed
	my $descs ="-";#Initialization
  
	# print "uniprotid: ".$uniprotId."\n";
  if ( (defined $res[0])  ){
    #print "pw_ids: ".$res[0]."\n"; #DEBUGCODE
    #my $ids = "";
    my @pwIds = split(",",$res[0]);#The pathways ids can be a list separated by comma
    #but these are database ids, we must find the corresponding names in the table
    foreach my $pwId (@pwIds){
      #print "Doing the query to extract seqId...\n";
      #print "Second Pathways query";
      my $select_handle2 = $dbh->prepare("SELECT description FROM pathwaysTable WHERE id='".$pwId."'");
      print "WARNING! In pathwaysInfo: Couldn't prepare second query; "
      unless defined $select_handle2; 
      $select_handle2->execute() or print $select_handle2->errstr;
      my @res2 = $select_handle2->fetchrow_array;
      
      #from this second result we can create the final string to add
      if (scalar(@res2) > 0){
        #This IF is created because sometimes something has already been
        #written. If this is the case then separate with comma
        #otherwise it is the first element
				
        if ($descs eq "-"){
          $descs = $res2[0];
        }else{
			$descs = $descs.$separator.$res2[0];
			
		}
      }
      #We need this else when the query gives an empty result. This means that there's no pathway  corresponding to the id
      else{
        print "ERROR in fetching informations from pathwaysTable in database. The id= $pwId did not correspond to anything.Check integrity of database!\n";
      }
    }
    #push ( @finalLine, $descs); 
		
  }#else {push ( @finalLine,"-");}
	
  push ( @finalLine, $descs); 
  return @finalLine;
}
=head2 pathwaysInfo

 Title : pathwaysInfo
 Usage : pathwaysInfo( -dbname => 'name of the GFF3 db',
			-use => 'username for dbname',
			-pass => 'password for dbname',
			-dbProt => 'the database uniprot created',
			-userProt => 'username for uniprot',
			-passProt => 'password for uniprot'
			 );

 Function: Returns the descriptions of the pathways.
			Takes in input an uniprot id that is usually the HSP name.
 
 Returns : an array

=cut
sub pathwaysInfoToCols{

  #VARIABLES FOR USE OF THE UNIPROT DATABASE
	my $uniprotId = shift;
	my $dbh = shift;
	my $viewPathway = shift;
  
  my @res = ();
  #Accessing database to extract enzyme ids corresponding to a certain uniprot id
	my $select_handle = $dbh->prepare("SELECT pw_ids FROM $viewPathway WHERE seqId='".$uniprotId."'");
	print "WARNING! In pathwaysInfo: Couldn't prepare first query; " unless defined $select_handle;
  $select_handle->execute()	or print $select_handle->errstr;
	@res = $select_handle->fetchrow_array;
							
	my @finalLine = ();#Initialization of the line that will be printed
  
	# print "uniprotid: ".$uniprotId."\n";
  if ( (defined $res[0])  ){
		my $pwL1 = ''; my $pwL2 = ''; my $pwL3 = '';
    #print "pw_ids: ".$res[0]."\n"; #DEBUGCODE
    #my $ids = "";
    my @pwIds = split(",",$res[0]);#The pathways ids is a list separated by comma
    #but these are database ids, we must find the corresponding names in the table
    foreach my $pwId (@pwIds){
      #print "Doing the query to extract seqId...\n";
      #print "Second Pathways query";
      my $select_handle2 = $dbh->prepare("SELECT description FROM pathwaysTable WHERE id='".$pwId."'");
      print "WARNING! In pathwaysInfo: Couldn't prepare second query; "
      unless defined $select_handle2; 
      $select_handle2->execute() or print $select_handle2->errstr;
      my @res2 = $select_handle2->fetchrow_array;
      
      #from this second result we can create the final string to add
      if (scalar(@res2) > 0){
        
        #The different pathways levels descriptions are separated by a ';' 
        #thus we split and construct three separate variables which will have the
        # default separator of Annocript to separate multiple values for same level
				my @pwLevels = split(";", $res2[0]);
				if (scalar(@pwLevels) > 0){
					$pwL1 .= $pwLevels[0].$separator;			
				}
				if (scalar(@pwLevels) > 1){
					$pwL2 .= $pwLevels[1].$separator;			
				}else{
					$pwL2 .= "-".$separator;#If a value
				}
				if (scalar(@pwLevels) > 2){
					$pwL3 .= $pwLevels[2].$separator;			
				}else{
					$pwL3 .= "-".$separator;
				}
      }
      #We need this else when the query gives an empty result. This means that there's no pathway  corresponding to the id
      else{
        print "ERROR in fetching informations from pathwaysTable in database. The id= $pwId did not correspond to anything. Check integrity of database!\n";
      }
    }
    
    #When the loop finishes, all the three variables are filled with the eventual multiple values
    #but we don't need the last separator
    if (scalar(@res) > 0){
		$pwL1 =~ s/(\Q$separator\E)$//;
		$pwL2 =~ s/(\Q$separator\E)$//;
		$pwL3 =~ s/(\Q$separator\E)$//;
	}   
	#All the strings are pushed inside the finalLine to return
	#if they exist
	if ($pwL1 ne "") {
      push(@finalLine, $pwL1);
    } else {push(@finalLine, "-");};
	if ($pwL2 ne "") {
      push(@finalLine, $pwL2);
    } else {push(@finalLine, "-");};
	if ($pwL3 ne "") {
      push(@finalLine, $pwL3);
    } else {push(@finalLine, "-");};
		
  }else {push ( @finalLine,"-","-","-");}#if the id has no pathway, three '-' are given
	
  return @finalLine;
}



=head2 separateBlastxTrSp

 Title : separateBlastxTrSp
 Usage : separateBlastxTrSp( -dbname => 'name of the GFF3 db',
			-use => 'username for dbname',
			-pass => 'password for dbname',
			-dbProt => 'the database uniprot created',
			-userProt => 'username for uniprot',
			-passProt => 'password for uniprot'
			 );

 Function: Takes features from an hash database and print them in a way that we choosed to represent:
          we print the TREMBL and SWISS-PROT result and the Organism and GO terms for the lowest evalue
						between the two.

 Returns : an array with the datas

=cut
sub separateBlastxTrSp{
	my $outHash= shift;
	my $id = shift;
	my $protDB = shift;
  my $NCflag = shift;
  
  #Variables to identify the type of database
	my $dbSP = 'sp';
	my $dbTR = 'tr';
	
	my @finalLine = ();#Initialization
   my @emptyLine = ();
   
   if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
	  #In empty line (no BLASTX/P results) the number of '-' depends from the number of informations related to BLASTX/P and actually are 24 using: 
	  #SP and TR results, the enzyme and Pathways columns. 
	  @emptyLine = ('-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-');	
   }
   else{
		#Here they are less since we are considering that the user did not ask for GO terms associated to proteins
	  @emptyLine = ('-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-');	
   }
	my $dbWinner = 'none';
	
	#There are 4 cases that we have to take in consideration.
	
	#The first case is when there are no results even with swiss-prot than with trembl
	if ( (!(defined $outHash->{$id}->{'sp'}) and !(defined $outHash->{$id}->{'tr'})) ){
			push (@finalLine,@emptyLine);
      #if there are no results the sequence could be non-coding
      #We increase twice: one for TR and another for SP
      $$NCflag++;#TR
      $$NCflag++;#SP
		#the second is when only sp has results
		}elsif (!(defined $outHash->{$id}->{'sp'}) and (defined $outHash->{$id}->{'tr'}) ){
				$dbWinner = 'tr';
        $$NCflag++;#TR
				#third, only trembl
				}elsif ((defined $outHash->{$id}->{'sp'}) and !(defined $outHash->{$id}->{'tr'})){
					$dbWinner = 'sp';
          $$NCflag++;#SP
				#but if each of them has a result.. then we evaluate the score
				#that better will describe the organism
				}elsif ($outHash->{$id}->{'tr'}->{'eval'} < $outHash->{$id}->{'sp'}->{'eval'}){
						$dbWinner = 'tr';
						}else{
						$dbWinner = 'sp';
						}
			
	my $partTime = time;
	my $partDuration = undef;
		 
	if ( $dbWinner ne 'none') {
    #We put in the resulting table only the results that have an evalue lower than a given threshold
		if ($outHash->{$id}->{$dbWinner}->{'eval'} <= $configHash->{'evalMax'} ){
			if ( $dbWinner eq 'sp'){
				my ($lineBackSp, $osSp) = blastxInfosUniprot($outHash,$id,$dbSP,'OS_YES', $protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfosSp occurs ".$partDuration." seconds\n";
				#$partTime = time;
			
				my ($lineBackTr, $osTr) = blastxInfosUniprot($outHash,$id,$dbTR,'OS_NO', $protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfosTr occurs ".$partDuration." seconds\n";
				#$partTime = time;
        
        #Printing the 
				push(@finalLine,@$lineBackSp);
				push(@finalLine,@$lineBackTr);
				if (defined($osSp)){
					push(@finalLine,$osSp);
				} else {
					push(@finalLine,'-');
				}
				if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
					push(@finalLine,blastxGOInfos($outHash,$id,$dbSP,$protDB));
			  }
				#$partDuration = time-$partTime;
				#print "blastxInfosGO occurs ".$partDuration." seconds\n";
				#$partTime = time;
			}elsif ( $dbWinner eq 'tr'){
				my ($lineBackSp, $osSp) = blastxInfosUniprot($outHash,$id,$dbSP,'OS_NO',$protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfosSp occurs ".$partDuration." seconds\n";
				#$partTime = time;
			
				my ($lineBackTr, $osTr) = blastxInfosUniprot($outHash,$id,$dbTR,'OS_YES',$protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfosTr occurs ".$partDuration." seconds\n";
				#$partTime = time;
					
				push(@finalLine,@$lineBackSp);
				push(@finalLine,@$lineBackTr);
					
				if (defined($osTr)){
					push(@finalLine,$osTr);
				}else{
					push(@finalLine,'-');
				}
				if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){			
					push(@finalLine,blastxGOInfos($outHash,$id,$dbTR,$protDB));
			  }
				#$partDuration = time-$partTime;
				#print "blastxInfosGO occurs ".$partDuration." seconds\n";
				#$partTime = time;
			}
		}else{push (@finalLine, @emptyLine);}	
   
	}
 	
	return @finalLine;
}

=head2 separateBlastxUfSp

 Title : separateBlastxUfSp
 Usage : separateBlastxUfSp( -outHash => 'hash containing the output from blast',
			-id => 'transcript id for which we are extracting results',
			-protDB => 'name of the protein database',
			-NCflag => 'flag for non coding information',
			-db2 => 'second db name (uniref or tr)',
			 );

 Function: Takes features from an hash database and print them in a way that we choosed to represent:
          we print the UniRef and SWISS-PROT result and the GO terms for the lowest evalue
						between the two.

 Returns : an array with the datas

=cut
sub printBlastProteinBestResults{
	my $outHash= shift;
	my $id = shift;
	my $protDB = shift;
  my $NCflag = shift;
  my $db2 = shift;
  
  
  #Variables to identify the type of database
	my $db1 = 'sp';
	#my $db2 = 'uniref';
	
	my @finalLine = ();#Initialization
   my @emptyLine = ();
   
   if ( ($configHash->{'goTermsAss'} eq 'proteins') or ($configHash->{'goTermsAss'} eq 'both')){
	  #In empty line (no BLASTX/P results) the number of '-' depends from the number of informations related to BLASTX/P and actually are 31 using: 
	  #SP and TR/uniref results, the enzyme and Pathways columns. 
	  @emptyLine = ('-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-');	
   }
   else{
		#Here they are less since we are considering that the user did not ask for GO terms associated to proteins (25)
	  @emptyLine = ('-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-');	
   }
	my $dbWinner = 'none';
	
	#Determine which result has the best e-value 
	#There are 4 cases that we have to take in consideration.
	#The first case is when there are no results even with swiss-prot than with trembl
	if ( (!(defined $outHash->{$id}->{$db1}) and !(defined $outHash->{$id}->{$db2})) ){
			push (@finalLine,@emptyLine);
      #if there are no results the sequence could be non-coding
      #We increase twice: one for TR/unrief and another for SP
      $$NCflag++;#TR/uniref
      $$NCflag++;#SP
		#the second is when only sp has results
		}elsif (!(defined $outHash->{$id}->{$db1}) and (defined $outHash->{$id}->{$db2}) ){
				$dbWinner = $db2;
        $$NCflag++;#TR/uniref
				#third, only trembl/uniref
				}elsif ((defined $outHash->{$id}->{$db1}) and !(defined $outHash->{$id}->{$db2})){
					$dbWinner = $db1;
          $$NCflag++;#SP
				#but if each of them has a result.. then we evaluate the score
				#that better will describe the organism
				}elsif ($outHash->{$id}->{$db2}->{'eval'} < $outHash->{$id}->{$db1}->{'eval'}){
						$dbWinner = $db2;
					}else{
						$dbWinner = $db1;
					}
			
	my $partTime = time;
	my $partDuration = undef;
		 
	if ( $dbWinner ne 'none') {
    #We put in the resulting table only the results that have an evalue lower than a given threshold
		if ($outHash->{$id}->{$dbWinner}->{'eval'} <= $configHash->{'evalMax'} ){
			  $os_tax = "os";
			  #print "Printing information for $id from $db1\n";#DEBUGCODE
			  my $lineBack1 = blastxInfosUniprot($outHash,$id,$db1,$protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfos1 occurred ".$partDuration." seconds\n";
				#$partTime = time;
				#print "Printing information for $id from $db2\n";#DEBUGCODE
				if ($db2 eq 'uniref'){$os_tax = "tax";}
				my $lineBack2 = blastxInfosUniprot($outHash,$id,$db2,$protDB);
				#$partDuration = time-$partTime;
				#print "blastxInfos2 occurred ".$partDuration." seconds\n";
				#$partTime = time;
				
			  #Printing the information about db1
				push(@finalLine,@$lineBack1);
				#Printing the information about db2
				push(@finalLine,@$lineBack2);
				#Printing GO info for the dbWinner
				push(@finalLine,blastxGOInfos($outHash,$id,$dbWinner,$protDB));
		}else{push (@finalLine, @emptyLine);}	
	}
	return @finalLine;
}

=head2 nonCodingInfo

 Title : nonCodingInfo
 Usage : nonCodingInfo( -outHash => 'the hash with data',
			-id => 'identifier of the transcript');

 Function: Takes features from an hash and print them in an array in the order to print on the final table.

 Returns : an array with the data

=cut
sub nonCodingInfo{
	#INPUT VARIABLES
	my $outHash = shift;
	my $id = shift;
  my $NCflag = shift;#flag to use to find out if the sequence is non-coding

	my $ncElement = "-";#Initialization
	
	#print "Extracting portrait informations...\n";
	if( defined $outHash->{$id}->{'codpot'} ){
			#push(@finalLine,$outHash->{$id}->{'codpot'});
      $ncElement = $outHash->{$id}->{'codpot'};
      if ( $outHash->{$id}->{'codpot'} <= $configHash->{'NCThresh'}){
        $$NCflag++;#if there are no results the sequence could be non-coding
      }
  }#else {push($ncElement,"-");}
  
  
	return  $ncElement;
}


=head2 orfInfo

 Title : orfInfo
 Usage : orfInfo( -outHash => 'the hash with data',
			-id => 'identifier of the transcript');

 Function: Takes features from an hash and print them in an array in the order to print on the final table.

 Returns : an array with the data

=cut
sub orfInfo{
	#INPUT VARIABLES
	my $outHash = shift;
	my $id = shift;
  my $NCflag = shift;#flag to use to find out if the sequence is non-coding

	my @finalLine = ();#Initialization
	
	#print "Extracting blastn informations...\n";
	if( defined $outHash->{$id}->{'dna2pep'}){	
			push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'l'});
			push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'os'});
			push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'of'});
      if ($outHash->{$id}->{'dna2pep'}->{'l'} <= $configHash->{'NCORFLength'}){
        $$NCflag++;#if there are no results the sequence could be non-coding
      }
  }else {push(@finalLine,"-","-","-");}
  
  
  
	return @finalLine;
}



=head2 countInfo

 Title : countInfo
 Usage : countInfo( -outHash => 'hash with data',
			-id => 'id of the feature',
			 );

 Function: Takes features from an hash and print them in an array in the order to print on the final table.

 Returns : an array with the data

=cut
sub countInfo{
	#INPUT VARIABLES
	my $outHash = shift;
	my $id = shift;
 			
	my @finalLine = ();#Initialization
	
	#print "Extracting blastn informations...\n";
	if( defined $outHash->{$id}->{'count'}){
			push(@finalLine,$outHash->{$id}->{'count'});
  }else {push(@finalLine,"-");}
	return @finalLine;
}


=head2 printORFFastaFile

 Title   : printORFFastaFile
 Usage   : printORFFastaFile(  -fastaOutFile => 'name of fasta file where write',
			-dbname => 'name of the database GFF',
			-user => 'username for the db',
			-pass => 'password for the db'
                         );

 Function: Prints on a fasta file the longest orf with strand and frame both for 
           the orf and for the trembl and swiss prot result.
          An example of header will be:
           >comp100008_c0_seq1 UniRef90_UPI00026531D8|1|1|DNA2PEP|177|+|1|
          >Transcript_1 SwissProt|O60488|1|2|Trembl|H2LVF5|1|2|DNA2PEP|721|+|2|

 Returns : nothing. It writes on a file in output in the working folder

=cut
sub printORFFastaFile{;

	my $fastaOutFile = shift; 
	my $outFolder = shift;
		    
  my $dbname = $configHash->{'gffDB'};
	my $user = $configHash->{'mySqlUser'};
  my $pass = $configHash->{'mySqlPass'};
    
	$fastaOutFile = $outFolder."/".$fastaOutFile;
	#Opening the gffdb connection
	my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
		                                 -user => $user,
		                                 -pass => $pass,
		                                 -dsn => $dsnGff,#"dbi:mysql:$dbname",
		                                 -create => 0);

	print "I'm writing now ORF informations in a file fasta...please wait...";
	#FOR EACH SEQUENCE RECORDED IN THE GFF DB
  my @ids = $db->seq_ids();
	#print "Sequence ids: \n";
	
  #Opening a file where print results
	open (OUT_FILE,">$fastaOutFile") || die ( "\nERROR [$?]: Cannot open fasta file to write: ?" );

	my $counter=1;	
		foreach my $id (@ids){
			print OUT_FILE ">$id $separator";
			
			#Print on a line informations retrieved from the blastx output
			print OUT_FILE join("|", blastxInfoORF($outHash,$id));
			print OUT_FILE "|";#A pipe will separate the data
			
			#Print on a line informations retrieved from the orf program output
			print OUT_FILE join("|", orfInfoORF($outHash,$id));
			#print OUT_FILE "\n";
			
			#print " ".$counter." ";
		}
	print "...DONE!\nThe output has been written on ".$fastaOutFile.".\n";
	close(OUT_FILE);
}


=head2 blastxInfoORF

 Title   : blastxInfoORF
 Usage   : blastxInfoORF(  -outHash => 'the hash with data',
			-id => 'identifier of the transcript');

 Function: Takes features from a hash and print them in a way that we choosed to represent.

 Returns : an array with the data

=cut
sub blastxInfoORF{
	my $outHash = shift;
	my $id = shift;
  
  #Variables to identify the type of database
	my $dbSP = 'sp';
	my $dbTR = 'tr';
  my $dbUniref = 'uniref';
	
	my @finalLine = ();#Initialization
  
  if( (exists $outHash->{$id}->{$dbUniref}->{'name'}) or
      (exists $outHash->{$id}->{$dbSP}->{'name'}) or 
      (exists $outHash->{$id}->{$dbTR}->{'name'}) ) {
    #Print out for the uniprot knowledgebase database
    if(exists $outHash->{$id}->{$dbSP}->{'name'}){
      push(@finalLine,'SwissProt');
      push(@finalLine,$outHash->{$id}->{$dbSP}->{'name'});
      push(@finalLine,$outHash->{$id}->{$dbSP}->{'strand'});
      push(@finalLine,$outHash->{$id}->{$dbSP}->{'frame'});
    }
    #Print out for the uniref database
    if(exists $outHash->{$id}->{$dbUniref}->{'name'}){ 
      #Here we don't need the UNIREF tag because it is included in the name
      push(@finalLine,split(/_/,$outHash->{$id}->{$dbUniref}->{'name'}));
      push(@finalLine,$outHash->{$id}->{$dbUniref}->{'strand'});
      push(@finalLine,$outHash->{$id}->{$dbUniref}->{'frame'});
    } 
    if(exists $outHash->{$id}->{$dbTR}->{'name'}){ 
      push(@finalLine,'Trembl');
      push(@finalLine,$outHash->{$id}->{$dbTR}->{'name'});
      push(@finalLine,$outHash->{$id}->{$dbTR}->{'strand'});
      push(@finalLine,$outHash->{$id}->{$dbTR}->{'frame'});
    }
  }else {push(@finalLine,"-","-","-");}
   
	return @finalLine;
}


=head2 orfInfoORF

 Title   : orfInfoORF
 Usage   : orfInfoORF (   -outHash => 'the hash with data',
                          -id => 'identifier of the transcript');

 Function: Takes features from an hash and print them in an array 

 Returns : an array with the data

=cut
sub orfInfoORF{
	#INPUT VARIABLES
  my $outHash = shift;
	my $id = shift;
        
  		
	my @finalLine = ();#Initialization
	
	#print "Extracting blastx informations for swiss prot and tremble...\n";
  if(defined $outHash->{$id}->{'dna2pep'}){
    push(@finalLine,'DNA2PEP');
		push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'l'});
		push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'os'});
		push(@finalLine,$outHash->{$id}->{'dna2pep'}->{'of'});
		#push(@finalLine,"\n");
			
		#EXPERIMENTAL
		if (defined $outHash->{$id}->{'dna2pep'}->{'orf'}){
      my $orf = $outHash->{$id}->{'dna2pep'}->{'orf'};
			my $seq="\n";
			while ($orf =~ /(.{1,60})/g) {
				$seq= $seq.$1."\n";
			}
			push(@finalLine,$seq);
		}else{push(@finalLine,'-');}
			#print "len:".$len[0]." strand: ".$strand[0]." frame: ".$frame[0]."\n";
  }else {push(@finalLine,"-","-","-","-");}
	return @finalLine;
}


1;
