#############################execute_programs.pm##############################
# This module is to execute programs and to extract the output of Annocript  #
##############################################################################

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

package PROGRAMS_EXEC::execute_programs;
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( execute_programs check_presence run_blast_indexing check_makeblastdb_out );

use strict;
use warnings;
use Data::Dumper;
use Cwd;#To change directories
 
#GENERAL UTILITY MODULES
use File::Copy;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use Term::ReadKey;
use Time::HiRes qw( time );
use IO::Handle;#To immediately print with autoflush 

#SEQUENCE MANAGEMENT MODULES
use Bio::SeqIO;
use Bio::SearchIO;#To use BioPerl sequence management
use Bio::AlignIO;

#DATA MANAGEMENT MODULES
use DBI;
use DBD::mysql;

#FOR PARALLEL EXECUTIONS
use Parallel::ForkManager;

#Annocript Libraries
use USEFUL::utilities qw(delete_file delete_directory exec_command 
                    count_sequences check_presence extract_name
                     extract_special_chars create_slices indexed_db_present
                     append_2_file);
                     

#EXTRACT PARAMETERS FROM THE CONFIG FILE
my $configHash;

#An array contains the arguments to be passed. If something that you need is not here please add	
my @blastArgs = ('word_size','evalue','num_threads','num_descriptions','num_alignments','max_target_seqs','outfmt','threshold','matrix');

#The maximum  and minimum lengths permitted for a sequence to be computed by portrait
my $maxPortraitSeqLength; 
my $minPortraitSeqLength;
my $maxPercNs;

#VARIABLE OF CONFIGURATION
my $fastaSeqs;
my $testFILE;
my $dataFolder;
my $outFolder;
my $gffFolder;
my $sessionFolder;
my $orfLEN ;
my $dbDataFolder;
  
#CODING/NONCODING CONFIG PARAMETERS /SEQUENCES WILL BE CLASSIFIED AS CODING OR NON-CODING AS PER THE PARAMETERS BELOW USED AS CUT-OFF)
my $NCSeqLength;
my $NCORFLength;
my $NCThresh;


##############################################READ CONFIGURATION PARAMETERS###############################################

#NEW INPUT VARIABLES
my $blastxPath;
my $blastnPath;
my $blastpPath;
my $tblastnPath;
my $rpstblastnPath;
my $rpsblastPath;
my $portraitPath;
my $dna2pepPath ;
my $makeblastdbPath;

#DATABASES PATHS
my $rfamDB;
my $swissProtDB;
my $tremblDB;
#my $uniprotKBDB;
my $cdDB;
my $cdDBPath;
my $query;
my $threads;
	
my $totSequences = 0;#Number of sequences in the file


      
=head2 execute_programs
 Title  : execute_programs
 Usage  : execute_programs(  );

 Function: 	This function will check and start each of the programs given in input an hash with parameters for their configuration
			and the query sequences file.
			
 Returns : nothing
=cut
sub execute_programs {
  $configHash = shift;
  
  STDOUT->autoflush(1);#This makes STDOUT hot in the sense that everything will be print immediately

  #Taking parameters in variables from the hash table		
  $testFILE = $configHash->{'fastaSeqs'};
	$sessionFolder = $configHash->{'sessionFolder'};
	$dataFolder = $sessionFolder."/".$configHash->{'dataFolder'};
	$outFolder = $sessionFolder."/".$configHash->{'outFolder'};
	$gffFolder = $sessionFolder."/".$configHash->{'gffFolder'};
	$fastaSeqs = $configHash->{'fastaSeqs'};
	$dbDataFolder = $configHash->{'dbDataFolder'};
	
  $orfLEN = $configHash->{'orfLength'};#10;
		
	#CODING/NONCODING CONFIG PARAMETERS /SEQUENCES WILL BE CLASSIFIED AS CODING OR NON-CODING AS PER THE PARAMETERS BELOW USED AS CUT-OFF)
	$NCSeqLength = $configHash->{'NCSeqLength'};#300;
	$NCORFLength = $configHash->{'NCORFLength'};#100;
	$NCThresh = $configHash->{'NCThresh'};#0.5;
  
  #Values to exclude sequences
  $minPortraitSeqLength = $configHash->{'min_port_seq_length'};
  $maxPortraitSeqLength = $configHash->{'max_port_seq_length'};
  $maxPercNs = $configHash->{'max_perc_ns'};
	##############################################READ CONFIGURATION PARAMETERS###############################################

	#NEW INPUT VARIABLES
	#These are the path to the programs that have to be present in the system
	$blastxPath = $configHash->{'blastxPath'};
	$blastnPath = $configHash->{'blastnPath'};
  $tblastnPath = $configHash->{'tblastnPath'};
  $blastpPath = $configHash->{'blastpPath'};
	$rpstblastnPath = $configHash->{'rpstblastnPath'};
  $rpsblastPath = $configHash->{'rpsblastPath'};
	$portraitPath = $configHash->{'portraitPath'};
	$dna2pepPath = $configHash->{'dna2pepPath'};
	$makeblastdbPath = $configHash->{'makeblastdbPath'};

	#DATABASES PATH
	$rfamDB = $configHash->{'ncDB'};
	$swissProtDB = $configHash->{'swissProtDB'};
	$tremblDB = $configHash->{'tremblDB'};
	#$uniprotKBDB = $configHash->{'uniprotKBDB'};
	$cdDB = $configHash->{'cdDB'};
  $cdDBPath = $configHash->{'cdDBPath'};
	$query = $configHash->{'fastaSeqs'};
	$threads = $configHash->{'threads'};
	
  #Variables for computing execution times
	my $globalStart = time;
	my $partTime  = time;
	my $partDuration = undef;
  
  my $programRuns = 0;
	################################################ Annotate Sequences #######################################################
	#Try if query file exists and opens
	my $seqsPath = $sessionFolder."/".$fastaSeqs;
	open(FH, $seqsPath) || die ( "\nERROR [$?]: No query file in $seqsPath!!!: ?" );
	close(FH);
	
	#$totSequences = `fgrep -c '>' $seqsPath`; #Extract number of sequences from the query file
	$totSequences = $configHash->{'totSequences'};
 
 	#Do Blastx
	if ( ($configHash->{'doBlastxSP'} eq "YES") or ($configHash->{'doBlastxTRorUf'} eq "YES")){
    my $blastToUse;
    my $blastToUsePath;
    
    #Here we check if the sequences are aminoacidic or peptidic because in the first case
    #blastx will be used, in the second one blastp
    if ($configHash->{'sequencesAlphabet'} eq 'dna'){
        $blastToUsePath = $blastxPath;
        $blastToUse = 'BLASTX';
      }else{
         $blastToUsePath = $blastpPath;
         $blastToUse = 'BLASTP';
        }
      
		if ( $configHash->{'dbInUse'} eq 'uniprotkb' ){
			if ($configHash->{'doBlastxSP'} eq "YES"){
				print "\n(".scalar(localtime)."): Running $blastToUse against SwissProt DB\n";#LOG
				my $blastxSPOutPath = 'blastx/'.$configHash->{'blastxSPOut'};
				&run_blast($fastaSeqs,$blastToUsePath,$configHash->{'swissProtDB'}, $threads,$blastxSPOutPath);
				#Check if the output is good
				$blastxSPOutPath = $blastxSPOutPath.'1';
				if (&checkOutput($blastxSPOutPath,$totSequences,"blast") == 0){
					print "WARNING! There had been some problems in the blastx running with swissprot. The number of sequences".
					"in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
				}
				$partDuration = time - $partTime;
				print "Finished $blastToUse against SwissProt in\t".$partDuration." seconds\n\n";
				$partTime = time;		
				$programRuns++;	
			}

      if ($configHash->{'doBlastxTRorUf'} eq "YES"){
				print "\n(".scalar(localtime)."): Running $blastToUse against Trembl DB\n";#LOG
				my $blastxTROutPath = 'blastx/'.$configHash->{'blastxTROut'};
				&run_blast($fastaSeqs,$blastToUsePath, $configHash->{'tremblDB'}, $threads,$blastxTROutPath);
				#Check the output if is good
				$blastxTROutPath = $blastxTROutPath.'1';
				if (&checkOutput($blastxTROutPath,$totSequences,"blast") == 0){
					print "WARNING! There had been some problems in the blastx running with trembl. The number of sequences".
					"in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
				}
				$partDuration = time - $partTime;
				print "Finished $blastToUse against TrEMBL in\t".$partDuration." seconds\n\n";
				$partTime = time;
				$programRuns++;
			}
    }
    if ( $configHash->{'dbInUse'} eq 'uniref' ){
			if ($configHash->{'doBlastxSP'} eq "YES"){
				print "\n(".scalar(localtime)."): Running $blastToUse against SwissProt DB\n";#LOG
				my $blastxSPOutPath = 'blastx/'.$configHash->{'blastxSPOut'};
				&run_blast($fastaSeqs,$blastToUsePath,$configHash->{'swissProtDB'}, $threads,$blastxSPOutPath);
				#Check if the output is good
				$blastxSPOutPath = $blastxSPOutPath.'1';
				if (&checkOutput($blastxSPOutPath,$totSequences,"blast") == 0){
					print "WARNING! There had been some problems in the blastx running with swissprot. The number of sequences".
					"in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
				}
				$partDuration = time - $partTime;
				print "Finished $blastToUse against SwissProt in\t".$partDuration." seconds\n\n";
				$partTime = time;
				$programRuns++;
			}
			if ($configHash->{'doBlastxTRorUf'} eq "YES"){
				print "\n(".scalar(localtime)."): Running $blastToUse against UniRef DB\n";#LOG
				my $blastxUnirefOutPath = 'blastx/'.$configHash->{'blastxUnirefOut'};
				&run_blast($fastaSeqs,$blastToUsePath,$configHash->{'unirefDB'}, $threads,$blastxUnirefOutPath);
				#Check if the output is good
				$blastxUnirefOutPath = $blastxUnirefOutPath.'1';
				if (&checkOutput($blastxUnirefOutPath,$totSequences,"blast") == 0){
					print "WARNING! There had been some problems in the blastx running with UNIREF. The number of sequences".
					"in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
				}
				$partDuration = time - $partTime;
				print "Finished $blastToUse against UniRef in\t".$partDuration." seconds\n\n";
				$partTime = time;
				$programRuns++;
			}
    }
  }
    
  #Do rps blast
	if ($configHash->{'doRpstblastn'} eq "YES"){
    my $blastToUsePath;
    my $blastToUse;
    my $cdDB = $configHash->{'cdDB'};
     
    #Here we check if the sequences are aminoacidic or peptidic because in the first case
    #blastx will be used, in the second one blastp
    if ($configHash->{'sequencesAlphabet'} eq 'dna'){
        $blastToUsePath = $rpstblastnPath;
        $blastToUse = 'RPSTBLASTN';
      }else{
         $blastToUsePath = $rpsblastPath;
         $blastToUse = 'RPSBLAST';
        }
		print "\n(".scalar(localtime).") Running $blastToUse against $cdDB db\n";#LOG
    my $rpstBlastnOutPath = 'rpstblastn/'.$configHash->{'rpstblastnOut'};
    
    my $threads = $configHash->{'threads4Parallel'};
    
    blast_parallel_run($fastaSeqs, $blastToUsePath, $cdDBPath, $threads, $rpstBlastnOutPath);
    
    $partDuration = time - $partTime;
    print "Finished $blastToUse against $cdDB db in\t".$partDuration." seconds\n\n";
    $partTime = time;
    $programRuns++;
	}
  
  if ($configHash->{'doBlastn'} eq "YES"){
      my $blastToUsePath;
      my $blastToUse;
      #Here we check if the sequences are aminoacidic or peptidic because in the first case
      #blastx will be used, in the second one blastp
      if ($configHash->{'sequencesAlphabet'} eq 'dna'){
        $blastToUsePath = $blastnPath;
        $blastToUse = 'BLASTN';
      }else{
         $blastToUsePath = $tblastnPath;
         $blastToUse = 'TBLASTN';
        }
         
    print "\n(".scalar(localtime).") Running $blastToUse against the non-coding RNA db\n";#LOG
    my $blastnOutPath = 'blastn/'.$configHash->{'blastnOut'};
    &run_blast($fastaSeqs,$blastToUsePath,$configHash->{'ncDB'}, $threads,$blastnOutPath);
    #Check the output if is good
    $blastnOutPath = $blastnOutPath.'1';
    if (&checkOutput($blastnOutPath,$totSequences,"blast") == 0){
      print "WARNING! There had been some problems in the $blastToUse running. The number of sequences".
           "in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
    }
    $partDuration = time - $partTime;
    print "Finished $blastToUse against ncRNAs ".$configHash->{'ncDB'}." in\t".$partDuration." seconds\n\n";
    $partTime = time;
    $programRuns++;
  }
  if ($configHash->{'sequencesAlphabet'} eq 'dna') {
    if ($configHash->{'doPortrait'} eq "YES"){
      print "\n".scalar(localtime).": Running PORTRAIT...\n";#LOG
      &portrait_split_and_run($fastaSeqs, $portraitPath,$configHash->{'portraitOut'});#Run portrait splitting the file in pieces
      #We don't check the output if is good at this point because it is done inside the function
      
      $partDuration = time - $partTime;
      print "Finished Portrait running in\t".$partDuration." seconds\n\n";
      $partTime = time;
      $programRuns++;
    }
    if ($configHash->{'doDna2Pep'} eq "YES"){
      print "\n".scalar(localtime).": Running dna2pep...\n";
      #Check the output if is good
      &run_d2p ($fastaSeqs, $dna2pepPath,$configHash->{'dna2pepOut'});
      if (&checkOutput($configHash->{'dna2pepOut'},$totSequences,"dna2pep") == 0){
        print "WARNING! There had been some problems in the dna2pep running. The number of sequences".
            "in the output differs from the number of sequences in the query. Please check and try to restart Annocript\n";
      }
      $partDuration = time - $partTime;
      print "Finished dna2pep running in\t".$partDuration." seconds\n\n";
      $partTime = time;
      $programRuns++;
    }
  }else{ print "Portrait and dna2pep will not be performed because you don't have DNA sequences...\n";}
  
  if ( $programRuns > 0 ) {
		$configHash->{'db_versions_file'} = $dataFolder."/".$configHash->{'db_versions_file'};
		#Remove if already present
		if ( -e $configHash->{'db_versions_file'} ){
				delete_file($configHash->{'db_versions_file'});
		}
		#Print useful information in the versions file
		my $info_string = "Annocript:".$configHash->{'AnnocriptVer'}."\n".
											"Annocript database:".$configHash->{'dbTablesVersion'}."\n".
											"UniProt version:".$configHash->{'dbVersion'}."\n".
											"UniProt database used:".$configHash->{'dbInUse'}."\n".
											"Domains database used:".$configHash->{'cdDB'}."\n";
		if ($configHash->{'dbInUse'} eq 'uniref'){
			$info_string .=	"UniRef type used:".extract_name($configHash->{'unirefDB'},1)."\n";
		}
											
		append_2_file($info_string,$configHash->{'db_versions_file'});
		#Calculate the time needed for all the processes
		my $globalDuration = time - $globalStart;
		print "\n".scalar(localtime).": Annocript has finished the execution of $programRuns programs in: $globalDuration seconds\n\n";#LOG
	}
}

 
############################CHECK OUTPUTS###########################

=head2 checkOutput
 Title  : checkOutput
 Usage  : checkOutput(  - fileOut => 'the file to check',
                        - numToCompare => 'the number to compare with the dimensions',
                        - algorithm => 'the algorithm used to have the output'
                        );

 Function: 	Check if the blast output have the same number of results as the number of queries in the fasta file.
			THIS CHECK IS POSSIBLE ONLY WHEN THE BLAST OUTPUT IS IN THE DEFAULT FORMAT
			OTHERWISE AN ERROR WILL BE REPORTED BUT THE CHECK IS WRONG
 Returns : 1 if check is ok, 0 otherwise

=cut
sub checkOutput{
	my $fileOut = shift;
  my $numToCompare = shift;
	my $algorithm = shift;
  
  
  my $ok = 0;#It's 1 when the numbers are equal
  
  #If the path has not already been built, builds it
  if ( !($fileOut =~ /$outFolder/) ) {
      $fileOut = $outFolder."/".$fileOut;#Add the path to the output folder
  }
	
	
	#THIS CHECK IS POSSIBLE ONLY WHEN THE BLAST OUTPUT IS IN THE DEFAULT FORMAT
	#OTHERWISE AN ERROR WILL BE REPORTED BUT THE CHECK IS WRONG
  if (-e $fileOut){
    if (!(defined $configHash->{'outfmt'}) or ( $configHash->{'outfmt'} eq "0")){
      if ($algorithm eq "blast"){
        my $seqInOut = `grep -c 'Query=' $fileOut`;#Count the number of queries done
        if ($seqInOut == $numToCompare){
          $ok = 1;
        }
      }
      if ($algorithm eq "portrait"){
        my $seqInOut = `grep -c '>' $fileOut`;#Count the number of queries done
        if ($seqInOut == $numToCompare){
          $ok = 1;
						}
				}
		if ($algorithm eq "dna2pep"){
			my $seqInOut = `grep -c '_rframe' $fileOut`;#Count the number of queries done
			if ($seqInOut == $numToCompare){
				$ok = 1;
			}
		}
    }else {$ok=1}
  }else{
    print "WARNING: File $fileOut doesn't exists...impossible to check consistency!\n";
  }
	return $ok;
}


#################################################### RUN PROGRAMS ##########################################################

=head2 run_blast_indexing

 Title  : run_blast_indexing
 Usage  : run_blast_indexing( -file => 'query.fasta',
                -program => 'name of the program'
			    );

 Function: 	Starts makeblastdb that indicizes the fasta database that blast applications have to use

 Returns : the blast db indicized will be saved in the working directory

=cut
sub run_blast_indexing {
  my $file = shift;
  my $program = shift; 
  my $programPath = shift;
  my $blastInfoProgram = shift;
  my $dbDataFolder = shift;
  
  #$file = $dataFolder."/".$file;
	 
  if(-z "$file"){
    die ("\nERROR [$?]: $file file is empty: ?\n");
  }
  if(!-e "$file"){
    die ("\nERROR [$?]: $file not exists: ?\n");
  }
  
  print "Running indexing... \n";
  my $extensions = '';
  my $outFILE = extract_name($file,0);
  
  #If the database is composed of peptides
  if(($program eq "blastx" ) or ($program eq "blastp") ){
    my $command = "$programPath -in $file -title $outFILE -dbtype 'prot'";
    print $command."\n"; 
    $extensions = 'blastIndexedExts';
    
    ( system ($command)  ) == 0
    or die "ERROR [$?]: an error occurred while running $programPath: $?";
  }
  
	#If the database is composed of nucleotides
	if( ($program eq "blastn") or ($program eq "tblastx")){ 
    my $command = "$programPath -in $file -title $outFILE -dbtype 'nucl'";
		print $command."\n"; 
		
    $extensions = 'ncRNAExts';
   ( system ($command)  ) == 0
    or die "ERROR [$?]: an error occurred while running $programPath: $?";
	}
	
	#Get the info on how many output files are created
	my $command = "$blastInfoProgram -db $file -info > ".$file.".info";
	#print $command."\n";
	( system ($command)  ) == 0
    or die "ERROR [$?]: an error occurred while running $blastInfoProgram: $?";
    	
	#Checking if all the files are ok ....
	if ( indexed_db_present(extract_name($file,0),$extensions,$dbDataFolder) == 0 ){
		die ("\nERROR [$?]: $programPath run not successful: ?\n");
	}
	
  #The following code is not tested
  #if($program eq "rpstblastn"){ 
		#my @elems = split("/",$file);
		#my $outFILE = $elems[scalar(@elems)-1];
		#my $command = "$programPath -in $file -title $outFILE -dbtype 'nucl'";
    #print $command."\n"; 
    	
    #( system ($command)  ) == 0
    #or die "ERROR [$?]: an error occurred while running $programPath: $?";
    
   #if( indexed_db_present(extract_name($file,0),"ncRNAExts",$dbDataFolder) == 0 ){
      #die ("\nERROR [$?]: $programPath run not successful: $?\n"); 
    #}
	#}
  print "....DONE!\n";
}

=head2 check_makeblastdb_out

 Title  : check_makeblastdb_out
 Usage  : check_makeblastdb_out( -filePath => 'the path of the indexed file');

 Function: 	Check if the indexing of a database has been well done. In particular is used the function blastdbcmd
          that, with the -info parameter, returns a series of informations about the indexing otherwise the string
            'BLAST Database error: No alias or index file found for nucleotide database [name.fasta] in search path [/home/francesco/test::]'

 Returns : the blast output will be saved in the working directory

=cut
sub check_makeblastdb_out {
  my $filePath = shift;
  my $programPath = shift;

  
  my $out = $filePath.".out";
  my $returnVal = 1;
  my $command = "$programPath -db $filePath -info -out $out";
  
  ( system ($command)  ) == 0
  or die "ERROR [$?]: an error occurred while running $programPath: $?";
  
  if ( -z $out){
    $returnVal = 0;
  }
  
  return $returnVal;
}

=head2 run_blast

 Title  : run_blast
 Usage  : run_blast( -file => 'the query file',
            -programPath => 'path to the blast program'
            -db => 'the db to use with blast'
            -threads => 'number of thread to use for blast'
            -outFile => 'the name of the output blast file'
			    );

 Function: 	This is a function to run a generic blast by using the parameters given using the function build_query_args

 Returns : the blast output will be saved in the working directory

=cut
sub run_blast {
  my $file = shift;
  my $programPath = shift;
  my $db = shift;
  my $threads = shift;
  my $outFile = shift;
 
  # $file is the query file that we decided to keep in the session folder
  $file = $sessionFolder."/".$file;
 
  #Split the path in elements from which we will take only the last that is the name
  my $program = extract_name($programPath,0);  
  
  $outFile = $outFolder."/".$outFile.'1';
  
  if(-z "$file"){ 
    die ("\nERROR [$?]: $program stops due to lack of sequences: ?\n");
  }
  if(!-e "$file"){  
    die ("\nERROR [$?]: $program stops due to lack of sequence file: ?\n");
  }
 
  print "Running ".$program."...";
  print "Out file: ".$outFile."\n";
   my $query = $programPath.' -query '.$file.' -db '.$db;
  if( ($program eq "blastx") or ($program eq "blastp") ){ 
    $query = $query.' '.build_query_args('X').' -out '.$outFile;  
  }
  elsif( ($program eq "rpstblastn") or ($program eq "rpsblast") ){ 
    $query = $query.' '.build_query_args('RPS').' -out '.$outFile;
  }
  elsif( ($program eq "blastn") or ($program eq "tblastn") ){ 
    $query = $query.' '.build_query_args('N').' -out '.$outFile;
  }
  
  print "The query is :".$query."\n\n";
  ( system ($query)  ) == 0
  or die "ERROR [$?]: an error occurred while running $programPath: $?";
    
  if(-z "$outFile"){
    print ("\nWARNING: $program output file is empty (maybe no results)\n");  
  }
  if(!-e "$outFile"){
    die ("\nERROR [$?]: $program run not successful: ?\n");
  }
  print "...DONE!! Results printed on $outFile\n";
}

	
=head2 build_query_args
 Title  : build_query_args
 Usage  : build_query_args( -programId => 'what is the program (X,N,RPS)' );

 Function: For all the BLAST programs we written only one subroutine: run_blast that simply executes 
            the given in input program by building a string with all the parameters inserted. 
            The parameters are written in the configuration file and when one of them is not present 
            (param1 = ) then it will not be added to the string.
           This is the procedure to build it: we used an array with all the names of the parameters that can be used. 
           It is coded in the starting of execute_program.pm 
           ('word_size','evalue','num_threads','num_descriptions','num_alignments','max_target_seqs','outfmt');
           in the configuration file these parameters are written with these name plus a letter that is linked with
            the program using (X,RPS,N). For example: word_sizeRPS, is the word_size for RPSTBLASTN program;
            
            Suppose that BLASTX is currently called for run. Then the subroutine build_query_args('X') is called
             with input parameter 'X'.

            In the subroutine an array will be used to find each of its elements. This is the array of names of parameters. ù
            The letter 'X' will be added to each of the names and the corresponding parameter will be searched in the configHash.

 Returns : The string built

=cut
sub build_query_args {
	my $queryArgs = '';
	my $programId=shift;

  foreach my $blastArg (@blastArgs){
    my $arg = $blastArg.$programId;
		if (defined($configHash->{$arg})){ 
			$queryArgs = $queryArgs.' -'.$blastArg.' '.$configHash->{$arg};
		}
	}
	return	$queryArgs;
}
	
	


=head2 portrait_split_and_run
 Title  : portrait_split_and_run
 Usage  : portrait_split_and_run( -file => 'the query file',
            -programPath => 'path to the blast program'
             -outFile => 'the name of the output blast file'
			    );
          
 Function: 	Executes an application that gives a probability of be non coding to each of the sequences.
            Since Portrait can interrupt stochastically its execution when the number of sequences 
            is very high or for other factors, than has been decided that if there are more sequences than an established
            threshold, all the query file is splitted in pieces to execute one at a time.
            

 Returns : the portrait output will be saved in the working directory

=cut
sub portrait_split_and_run {
  my $startfile = shift;
  my $programPath = shift;
  my $outRealName = shift;
  
  #Erase special chars	from the principal fasta file and takes a new query file
  my $file = portrait_cleanFasta($startfile);
 
  #print "file: $file - startfile: $startfile - sessionFolder:  $sessionFolder\n";#DEBUGCODE
  # $file is the query file that we decided to keep in the session folder
  $file = $sessionFolder."/".$file;
  my $outFile = $file."_results_all.scores";
  
  #Split the path in elements from which we will take only the last, that is the name
  my @pathElems = split("/",$programPath);
  my $program = pop(@pathElems);

  print "Running ".$program."...\n";
  if(!-e "$file"){
    print ("\nWARNING: $program stops due to lack of sequence file ($file)\n");
  }else{
    if(-z "$file"){
      print ("\nWARNING: $program stops due to lack of sequences ($file)\n");
    }else{
      my $binDir = join("/",@pathElems);
      my $workingDir = getcwd;  

      #This is the name for each of the portrait outputs 
      my $mergedOut = $outFolder."/".$outRealName;
     
      #Delete any precedent version of output of Portrait
      #print "Deleting any precedent version of PORTRAIT output\n";#DEBUCODE
      delete_file($mergedOut);
      
      chdir $binDir;#Enters in the directory of PORTRAIT
      my $toExam = $configHash->{'split_num'};# sequences per time to analyze
      if ($totSequences > $toExam){ 
        
        #Here we put a control on the execution 
        print "..splitting the file in pieces with $toExam sequences...\n";
        portrait_exec_blocks($file,$programPath,$mergedOut,$toExam,$binDir);
        if(!(-e "$mergedOut")){
          # unlink($file."_results.orf",$file."_results_ncRNA.fas",$file.".log");
          print ("\nWARNING: $program run not successful! Can't find file $mergedOut. Annocript will continue without...\n");
        } 
        chdir $workingDir;
      }else{
        my $command = "$programPath -i $file -q 1";
        print "\n The command is: $command\n";
       
        ( system ($command)  ) == 0
        or die "ERROR [$?]: an error occurred while running $programPath: $?";
    
        chdir $workingDir; 
        if(!(-e "$outFile")){
          unlink($file."_results.orf",$file."_results_ncRNA.fas",$file.".log");
          print ("\nWARNING: $program run not successful! Can't find file $outFile. Annocript will continue without...\n");
        }
        rename $file."_results_all.scores",$outRealName;#changing the name of the portrait output
        move($outRealName,$outFolder) or print "ERROR [$!]: unable to move $outRealName in $outFolder\n";#Moving the output of Portrait to the output folder
        unlink($file."_results.orf",$file."_results_ncRNA.fas",$file.".log");
      }
      #unlink($file."_results.orf",$file."_results_ncRNA.fas",$file.".log");#DEBUGCODE
      #print "...DONE!\n";
    }
  }
  
  #if the user doesn't want to keep the file generated for Portrait
  #per default it will be removed
  if ( $configHash->{'keepTempFiles'} eq "NO"){
	  delete_file($file); 
  } 
}


=head2 portrait_split_and_run
 Title  : portrait_split_and_run
 Usage  : portrait_split_and_run( -seqFile => 'the query file',
                                  -programPath => 'path to the blast program'
                                   -mergedOut => 'the name of the output portrait file with all the pieces of output',
                                   -toExam => 'number of seqeucence in the block',
                                   -binDir => 'directory of portrait program'
                                );

 Function: 	THIS EXECUTES PORTRAIT SPLITTING THE FILE IN DIFFERENT PIECES THEN IT WILL RUN EACH 
              OF THEM SEPARATELY IN THIS WAY PORTRAIT WILL NOT HAVE STOCHASTICAL INTERRUPTIONS

 Returns : the portrait output will be saved in the working directory

=cut
  

sub portrait_exec_blocks {
  my $seqFile = shift;
  my $programPath = shift;
  my $mergedOut =shift;
  my $toExam = shift;
  my $binDir = shift;
    
  
  open(IN, $seqFile) or 
    die "ERROR: some problem occurred while opening $seqFile. Check permissions to access the folder and restart Annocript.\n";
  
  my $piecesFile = "sequence.fasta";
  open(SEQ_FILE,">$piecesFile") or 
    die "ERROR: some problem occurred while opening $piecesFile. Check permissions to access the folder and restart Annocript.\n";

  #Variables to control the blocks
  my $seqCount= 0;#Counts the sequences for the block
  my $stopSeq = 0;#    
  my $done = 0;#Says when the execution is succesfull, if faults are too much it will be 1 the same 
  
  #Variables to control the effective functioning of Portrait
  my $badRuns = 0;#Take in count the number of stochastic faults of Portrait
  my $maxBadRuns = 10;#Number of times that Annocript should restart portrait  for its stochastic faults
  my $stopPortrait = 0;#This is a flag to control and see if Portrait has to stop completely
  
  #my $totSequences = `fgrep '>' $seqFile | wc -l `; #Extract number of sequences from the query file
	
  #The execution of portrait is divided in blocks of $toExam sequences
  #the last block comprise all the remaining sequences
  #but if the total number of sequences is minor than $toExam, no division will be done
    
  my $blockNum = 0;
  
  my $command = "$programPath -i $piecesFile -q 1";
  print "\n The command is: $command\n";
  
  ( system ("touch $mergedOut") ) == 0
        or print "ERROR in execute_programs::portrait_exec_blocks: an error occurred while creating $mergedOut: $?";
        
  #We read all the file with the sequences
  while( (my $row = <IN>)  ) {
    if ($stopPortrait == 0){
        #If the line starts with a > then it will save only the string until the first space
      if($row =~ /^\>/){
        $stopSeq++;
        if ($toExam == 1) {print "Elaborating sequence: $stopSeq\n";}#DEBUGCODE
        #If the number of sequences to examine is reached and the block is not the last one
        if (($seqCount ==  $toExam) and ($blockNum < (int ( $totSequences/$toExam)) ) ){
         if ($toExam == 1) { print "sequenza: $stopSeq blocchi da:$toExam blocco num:$blockNum non superare:".int($totSequences/$toExam)."\n";}#DEBUGCODE
          #Close the file with the just wrote sequences
          close(SEQ_FILE);
          $blockNum++;#Give a number to the block
          #Execute the run until the output is obtained
          while ( $done == 0){
            ( system ($command)  ) == 0
            or print "ERROR in execute_programs::portrait_exec_blocks: an error occurred while running $command: $?\n";
          
            my $outFile = "sequence.fasta_results_all.scores";
            my $resultNum = `fgrep -c '>' $outFile`;#Extract the number of sequences outputed by portrait
            
            ##Copying the file in a fasta file with the number of block
            #( system('cp sequence.fasta sequence_'.$blockNum.'.fasta')  ) == 0
            #or print "ERROR in execute_programs::portrait_exec_blocks: an error occurred while copying $piecesFile: $?";
            
            #If the output has been produced (sequence.fasta_results_all.scores)
            if( (-e "$outFile") ){
              #...and the number of results is at least half of the sequences in the block
              if ( $resultNum >= $toExam/2 ){
                $done = 1;#This running is done
                $badRuns = 0;
                #Copy the output in mergedOut appending it
                ( system ("cat $outFile >> $mergedOut")  ) == 0
                or die "ERROR [$?]: an error occurred while appending $outFile to $mergedOut: $?\n";
                if ($toExam == 1) {print "Good! Number of elements in $blockNum output :".`fgrep -c '>' $outFile`."\n";}#DEBUG CODE
              }else{   
                $badRuns++;
                if ($badRuns >= $maxBadRuns){
                   $stopPortrait = 1;
                   $done = 1;
                   print ("\n[>>> FATAL ERROR: PORTRAIT run not successful! Annocript will continue without its ouput.<<<]\n");
                }else{ 
                  print ("\n[>>> WARNING: PORTRAIT is resolving very few sequences only. \nRestarting ($badRuns trial)....<<<]\n"); 
                }
              }
            }else{
              $badRuns++;
              if ($badRuns >= $maxBadRuns){
                $stopPortrait = 1;
                $done = 1;
                print ("\n[>>> FATAL ERROR: PORTRAIT run not successful! Annocript will continue without its ouput.<<<]\n");
              }else{ 
                print ("\n[>>> WARNING: PORTRAIT block $blockNum run not successful! Can't find file $outFile. Restarting ($badRuns trial)....<<<]\n");
              }
            }
            delete_file($binDir.'/sequence.fasta_*');
          }
          open(SEQ_FILE,">$piecesFile");
          print SEQ_FILE $row;
          $seqCount = 0;
          $done = 0;
        }else{
          print SEQ_FILE $row;
        }
        $seqCount++;
      }else{
        print SEQ_FILE $row;
      }
    }
  } 
  close(SEQ_FILE);  
  close(IN);  
  delete_file($binDir.'/sequence.fasta_*'); 
  
  my $remainingSeqs = `fgrep -c '>' $piecesFile`;
  if ( ($remainingSeqs > 0) and ($stopPortrait == 0) ){
    print "Completed the print the output for the blocks! \n Starting with the final piece of $remainingSeqs seqs...\n";
    #The last group of genes will not be 100 but they have to be
    #analyzed the same
    $done = 0;
        
    while ( $done == 0){
       ( system ($command)  ) == 0
                or print "ERROR [$?]: an error occurred while using $programPath: ?";
      my $outFile = "sequence.fasta_results_all.scores";
      if(-e "$outFile"){
        $done = 1;
        ( system ("cat $outFile >> $mergedOut")  ) == 0
                or print "ERROR [$?]: an error occurred while concatenating  $outFile and mergedOut: ?";
      }else{
        print ("\n[>>> WARNING: PORTRAIT block $blockNum run not successful! Can't find file $outFile. Restarting....<<<]\n");
      }
      delete_file($binDir.'/sequence.fasta_*');
    } 
    delete_file($binDir.'/sequence.fasta_*');
    print ("removed output files...\n");
  }else{print "Last block with zero sequences...\n";}
  
  return $mergedOut;
}
  
  
  
    
=head2 run_d2p
 Title  : run_d2p
 Usage  : run_d2p( -file => 'the query file',
            -programPath => 'path to the blast program'
           -outFile => 'the name of the output blast file'
			    );

 Function: 	Executes an application that gives a probability of be non coding to each of the sequences

 Returns : the dna2protein output will be saved in the working directory

=cut
sub run_d2p {
  my $file = shift;
  my $programPath = shift;
  my $outFile = shift;

  my $mode = $configHash->{'d2pMode'};
  #Since dna2pep has problems with headers. We have to completely delete all the informations
  #and retain only the name
  my $cleanFile = cleanHeader( $file);
 
  move($cleanFile,$sessionFolder) or print "ERROR [$!]: unable to move $cleanFile in $sessionFolder. Please check permissions.\n";
  # $file is the query file that we decided to keep in the session folder
  $cleanFile = $sessionFolder."/".$cleanFile;
  $outFile = $outFolder."/".$outFile;

 
  #Split the path in elements from which we will take only the last that is the name
  my @pathElems = split("/",$programPath);
  my $program = $pathElems[@pathElems-1];
 
  if(-z "$cleanFile"){
    print ("\nERROR [$?]: $program stops due to lack of sequences: ?\n");
  }
  if(!-e "$cleanFile"){
    print ("\nERROR [$?]: $program stops due to lack of sequence file: ?\n");
  }
  #Here we run DNA2PEP where we specify to obtain the tabular output (-O tab)
  # -o none - is to not use the start codon
  # -r all is to specify to read all the frames 
  my $command = "$programPath -O tab -o $mode -r all $cleanFile > $outFile";
  print "The command is: $command\ns";
  ( system ($command)  ) == 0
     or print "ERROR: system error occurred while running $programPath: $?";
 
  if(-z "$outFile"){
    print ("\nERROR [$?]: $program run not successful: ?\n");
  }
  if(!-e "$outFile"){
    print ("\nERROR [$?]: $program run not successful: ?\n");
  }
  delete_file($cleanFile);
  print "...DONE!\nDNA2PEP results written on $outFile\n";
}

##############################################UTILITIES#############################################


=head2 portrait_cleanFasta

 Title  : portrait_cleanFasta
 Usage  : portrait_cleanFasta( );

 Function: checks if in the input query file are there special chars.
			If yes, a code will erase them. This is done for the portrait execution.
			Then cleans the header to take only the compound name
			
 Algorithm: 
			- check for special chars from the query file with a command
			- if the special_char.txt is not empty I catch all the special chars 
			- with all the chars found I'll build a string that I'll put in a regex
			 Each of the special chars will be substituted by a "N".			
			
			
			- if the file has some other informations other than the name, in the header
			  then it will save only the string until the first space
			
      - Erase all the longest than 100 thousands sequences
        or portrait will give "Illegal division by zero"
 Returns : nothing

=cut
sub portrait_cleanFasta{
  my $queryFile = shift;
	
  my $portraitName = "portrait_".$queryFile;
	
	#Recognize special characters and put them inside a file
	my $spCharsFile = "special_char.txt";
	my $queryFile2;
    
  #Clean the file for protrait erasing all the useless informations from the sequences names
  my $cleanHQueryFile = cleanHeader($queryFile);	
  print "PORTRAIT CLEAN: Cleaned the header of the input to be accepted by PORTRAIT!\n";
  
  my $queryFilePathCompl = $cleanHQueryFile;#$sessionFolder."/".$cleanHQueryFile;
   
  
  #fgrep with -v inverts the sense of matching. So, for the lines that hasn't a '>'
	#sed simply apply the transformation of the regexp to all the file
	# The result is a file with only letters different from 'ATGCN'

  # ( system ("fgrep -v '>' $queryFilePathCompl|sed 's/[ATGCN]//ig'  > $spCharsFile")  ) == 0 # OLD COMMAND
  ( system ("fgrep -v '>' $queryFilePathCompl|sed 's/[ATGCN]//ig' | tr -d '\r \n' > $spCharsFile")  ) == 0
              or die "ERROR [$?]: an error occurred while extracting the special chars. (Using system commands: fgrep and sed: $?";
  print "PORTRAIT CLEAN: Checked if are present strange chars!\n";
              
	open(SP,"<$spCharsFile");
	#Extracting an array with all special chars
	my $spChars = "";
	my @spCharsArray = ();
	while (<SP>){
    my @splitted = split //;
    foreach my $spChar (@splitted){
      #If the char is not a space, a newline or is already present in the array
      if (($spChar ne " ") and ($spChar ne "\n") and ( (grep(/\Q$spChar\E/,@spCharsArray)==0) ) ){
        #print "Inserting $spChar\n";
        push(@spCharsArray,$spChar);
      }
    }
  }
  close(SP);	
  #Cancel the file used to memorize special chars 
  unlink($spCharsFile);
		
	#If special chars are present then create a string and match it with the clean file	
	if (@spCharsArray>0){
    foreach my $spChar (@spCharsArray)	{
      #print "Found: $spChar";
      $spChars = $spChars.$spChar;
    }
    $spChars =~ s/\s//g;
    print "\tUnknown chars found: ".$spChars."! They will be substituted with N. Only for PORTRAIT execution.\n";
        
    #Subsituting all the special chars with N and creating a different file for portrait
    open(PORT_QUERY,">".$portraitName);
    open(QUERY,"<$queryFilePathCompl");
        
    while(<QUERY>){
      my $line = $_;
      chomp($line);#Erase the newline char 
      if($line =~ m/^>/){
        print PORT_QUERY $line."\n";
      }else{
        $line =~ s/[$spChars]/N/g;
        print PORT_QUERY $line."\n";
        #print "$line\n";
      }
    }
    close(PORT_QUERY);
    close(QUERY); 
    $queryFile2 = $portraitName ;
	}else{
    print "\tNo special chars found in the query file!\n";
    $queryFile2 = $cleanHQueryFile;
  }
	
	  
  #Erase all the longest than 100 thousands sequences
  #or portrait will give "Illegal division by zero"
  print "PORTRAIT CLEAN: Checked the sequences! (Longest than $maxPortraitSeqLength are not accepted)\n";
  my $modQueryFile = portrait_correct_seqs($queryFile2);#$portraitName);
  rename $modQueryFile,$portraitName;
		
  #The last output has to stay in the session folder to be executed
  move($portraitName,$sessionFolder) or print "ERROR [$!]: unable to move $portraitName in $sessionFolder. Please check permissions.\n";
    
		
  print "PORTRAIT CLEAN: completed adaptation of your query file for PORTRAIT!\n";
  #Portrait will have another file for sequences
  return $portraitName;		
}
	
 
  
  
=head2 cleanHeader

 Title  : cleanHeader
 Usage  : cleanHeader( );

 Function: this function erase from the header of the sequences the informations that makes the name to be too long
			
 Returns : the name of the cleaned file

=cut		
sub cleanHeader {
  my $queryFile = shift;
	
	my $cleanFastaFile = "cleanh_".$queryFile;
  $queryFile = $sessionFolder."/".$queryFile;
	
	#CLEAN HEADER
  #When a fasta file has the name of the sequence plus other informations separated by some charachters
  #this script will give us a file cleaned of these
  open(OUT,">$cleanFastaFile ");
  open(IN,$queryFile);#OPEN PARSE FILE
 
  while(my $row = <IN>) {
    #If the line starts with a > then it will save only the string until the first space
    if($row =~ /^\>/){
      my @pieces = split(" ",$row);
      print OUT $pieces[0]."\n";
    }else{
      print OUT $row;
    } 
  }
  close(IN);
  close(OUT);

	return $cleanFastaFile;
}	



=head2 portrait_correct_seqs

 Title  : portrait_correct_seqs
 Usage  : portrait_correct_seqs( );

 Function: - this function erase the longest then $maxPortraitSeqLength sequences from the file
            - erase the shortest than $minPortraitSeqLength
            - erase the sequences with more then 20% of Ns
			
 Returns : nothing

=cut	
sub portrait_correct_seqs{
	my $queryFile = shift;
	
	my $output= "cutted_".$queryFile;
	#$queryFile = $sessionFolder."/".$queryFile;
		
	my $in = Bio::SeqIO->new(-format => 'fasta', -file=> $queryFile);
	my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$output");
           
	my $countLongest = 0;	
  my $countShort = 0;  
  my $countNs = 0;   
  my $remainingSeqs = 0;  
  
  while ( my $seq = $in->next_seq()) {
    my $tr_name = $seq->id;
    my $description = $seq->desc();
    my $str = $seq->seq; # get the sequence as a string
    my $basecount = $seq->length;
   
    #Counting the Ns using the () list context
    my $N = 'N';   
    my $Ns = () = $str =~ /$N/ig;  # $Ns is the number of N inside the sequence
   
    #Maximum admitted number of Ns in a sequence is maxPercNs of its length
    if ($Ns < ( ($seq->length/100) * $maxPercNs) ){
      if ($basecount > $minPortraitSeqLength) {
       if ($basecount < $maxPortraitSeqLength){
          my $seqrepeat = Bio::Seq->new(-seq => $str,
                        -display_id => $tr_name,
                        -description =>$ description); 
          $out->write_seq($seqrepeat);
          $remainingSeqs++;   
       }else{
          $countLongest++;
          print "Sequence $tr_name erased because it's longer ($basecount) than the maximum ($maxPortraitSeqLength)\n";
       }
      }else{
          $countShort++;
          print "Sequence $tr_name erased because it's shorter ($basecount) than the minimum ($minPortraitSeqLength)\n";
       }
    }else{
      $countNs++;
      print "Sequence $tr_name erased because it has $Ns Ns and it's length is $basecount (too much Ns!!!)\n";
    }
	}
  if (  ($countLongest >0) or ($countNs >0) or ($countShort >0) ){ 
    print "Found $countLongest sequences longest then $maxPortraitSeqLength nucleotides, $countShort shorter than $minPortraitSeqLength, ".
        "$countNs with more then".
    " 20% of Ns. They have been erased only for PORTRAIT. Good sequences are: $remainingSeqs\n".
    " No other programs will be influenced by this change.\n";
  }
  return $output;
}
##########################################PARALLEL EXECUTION##########################


=head2 blast_parallel_run

 Title  : blast_parallel_run
 Usage  : blast_parallel_run( - fasta => 'the name of the fasta file',
                              - programPath => 'the path to the program',
                              - db => 'the path of the db used by the program',
                              - processes => 'number of processes to be used',
                              - outFile => 'name of the output file' 
                              );

 Function: builds the command to be run, calls a routine to split the fasta file and then a routine for
          the execution of the processes.
			
 Returns : nothing

=cut	
sub blast_parallel_run{
  my $fasta = shift;
  my $programPath = shift;
  my $db = shift;
  my $processes = shift;
  my $outFile = shift;
  
   # $file is the query file that we decided to keep in the session folder
  my $sliceName = $fasta;#This will be the initial part of each slice
  $fasta = $sessionFolder."/".$fasta;
  $outFile = $outFolder."/".$outFile;
  
  my $workingDir = getcwd;
  
  #print "Working directory: $workingDir\n";#DEBUGCODE
  #Split the path in elements from which we will take only the last that is the name
  my @pathElems = split("/",$programPath);
  my $program = $pathElems[@pathElems-1]; 
  
  
  my $cmd = $programPath.' -db '.$db;
  #This works actually only for rps because blastn and blastx can be started with num_threads
  if( ($program eq "rpstblastn") or ($program eq "rpsblast") ){ 
    $cmd = $cmd.' '.build_query_args('RPS').' ';
  }
  #elsif( ($program eq "blastx") or ($program eq "blastp") ){ 
    #$cmd = $cmd.' '.build_query_args('X').' ';  
  #}
  #elsif( ($program eq "blastn") or ($program eq "tblastn") ){ 
    #$cmd = $cmd.' '.build_query_args('N').' ';
  #}
  
  #Build a name with 'slice' in the end
  $sliceName =~ s/\..+$//;#erase the .fasta extension from the name
  $sliceName .= '_slice';#add _slice
  
  #The slices are created from the $fasta file in a $sliceName folder
  create_slices($fasta,$sliceName,$processes, $sessionFolder);
  
  my $actualDir = getcwd;
  
  #Executes the command $cmd for all the slices in the folder
  parallel_execution($sliceName,$processes,$cmd,$outFile);
  
  chdir $workingDir;
  
  #if the user doesn't want to keep the temporary files, the folder with the slices are removed  
  #For default all the temporary files are removed
  if ( $configHash->{'keepTempFiles'} eq "NO"){
	  delete_directory($sessionFolder.'/'.$sliceName); 
  }
}
  
  
=head2 parallel_execution

 Title  : parallel_execution
 Usage  : parallel_execution( - fastaFolder => 'folder with the fasta files',
                              - processes => 'number of processes to be used',
                              - cmd => 'the command to be used. It has to be without query and out',
                              - out => 'path of the output file' 
                              );

 Function: for each FASTA file present in the fastaFolder executes once the command cmd and creates a different output
          whose name is given by $out and a progressive number.
			
 Returns : nothing

=cut	 
sub parallel_execution{
  my $fastaFolder = shift;
  my $processes = shift;
  my $cmd = shift;#This will be only the string with parameters, no in and out files
  my $out = shift;#This comprises the rps blast out path
  
  print "Fetching files from $fastaFolder and executing:\n $cmd\n";# I'll write a log in $log.\n"; 
  #Creates an object to fork processes in parallel
  my $pm = new Parallel::ForkManager($processes);

  #You can define a subroutine which is called when a child is terminated. It is called in the parent process.
  #  - pid of the process, which is terminated
  # - exit code of the program
  $pm->run_on_finish(
    sub {
      my($pid,$exit_code) = @_;
      #print "** Just got out of the pool with PID $pid and exit code: $exit_code\n";#DEBUCODE
    }
  );

  #print "Working directory inside parallel_execution: $actualDir\n";#DEBUCODE
  
  #takes all the fasta files from the folder
  my @fasta = glob('*.fasta');
  #print "Folder just open. We found: ".scalar(@fasta)." files...\n";#DEBUCODE
  my $numFasta = 1;
  #Starts the processes and prints the PID of the childs
  foreach my $fasta (@fasta) {
    
    my $outName = $out.$numFasta;
    #print "Working on File: $outName\n";#DEBUCODE
    #Build the command to run
    my $run = $cmd.'-query '.$fasta.' -out '.$outName;
    sleep 1;
    $numFasta++;
    # Forks and returns the pid for the child:
    my $pid = $pm->start and next;
   
    # Here is the parallelized block
    # -----------
    #print "$pid ---> running on $fasta\n";
    #print  "Launching $run\n";
    print "-query ".$fasta." -out ".$outName."\n";
    exec_command($run);
    
    # We check if the output is ok
    my $numToComp = count_sequences($fasta);
    if (&checkOutput($outName,$numToComp,"blast") == 0){
			print "WARNING! There had been some problems in the running. The number of sequences in $outName".
					"differs from the number of sequences in the query. Please check and try to restart Annocript\n";
    }
    # Terminates the child process
    $pm->finish;
    print " ...DONE! Results printed on $out(i)\n";
  }
  $pm->wait_all_children; 
}

1;
