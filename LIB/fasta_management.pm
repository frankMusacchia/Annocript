
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
    
package LIB::fasta_management;
## fasta_management.pm
#Author: Francesco Musacchia  (2016)
#Permits the management of an files with many manipulation functions
BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( checkFastaFormat detect_fasta_type count_sequences create_slices);
}
use strict;
use warnings;
#SEQUENCE MANAGEMENT MODULES
use Bio::SeqIO;
use File::Copy;#To manage files



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
  
  die "Unable to change directory in $slicesDir.\n" unless chdir $slicesDir;
  
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

=head2 checkFastaFormat

 Title   : checkFastaFormat
 Usage   : checkFastaFormat( - fastaFilePath = path of the fasta file
                             - suggestedLength = suggested length for the line
                               );

 Function:  this subroutine controls the fasta file. It should respect
				the following rules:
				 - all the lines, header included, should be shorter than 80 chars
				 - no spaces between the symbol '>' and the transcript name
				 - IUPAC characters are only allowed (http://www.bioinformatics.org/sms/iupac.html)
        To check for the IUPAC code will be detected before if sequences
        are nucleotides or aminoacids.
        Line numbers in the file with the produced problem are given
            
 Returns :  Number of warnings obtained

=cut
sub checkFastaFormat{  
  my $fastaFilePath = shift;
	my $suggestedLength = shift;
	my $nuclIUPAC = shift;#"ACGTURYSWKMBDHVN\.\-";
	my $protIUPAC = shift;#"ACDEFGHIKLMNPQRSTVWY";
		
	#Variables for the subroutine
	my $totSequences = 0;
	my $lineNum = 0;
	my $warnings = 0;
	
	#Here I use variables to detect problems in the lines and their positions
	#in such a way that all the position are concatenated 
	my $tooLongHeaders = 0;#Number of headers too long
	my $tooLongHPos = '';
	my $tooLongSeqs = 0;#Number of sequences too long
	my $tooLongSPos = '';
	my $spaceNearName = 0;#Number of sequences wiht spaces after >
	my $spaceNearNamePos = '';
	my $noIUPACRows = 0;#Number of rows with non IUPAC characters
	my $noIUPACRowsPos = '';
		
	#Detecting if in the Fasta there are nucleotides or peptides
	my $seqType = detect_fasta_type($fastaFilePath);
	my $IUPACCode;
	
	#print "Detected sequence type: $seqType\n";#DEBUGCODE
	if ( ($seqType eq 'dna') or ($seqType eq 'rna') ){
		$IUPACCode = $nuclIUPAC;
	}else{
		$IUPACCode = $protIUPAC;
	}
	
	#Open the fasta
	open (FASTA,"<$fastaFilePath");
	
	#Go through the file to look for problems..
	while (my $row = <FASTA>){
			$lineNum++;
			#Here I check the header length
			if ($row =~ /^>/){
					chomp($row);
					$totSequences++;
					if (length($row) > $suggestedLength){
						$tooLongHeaders++;
						if ($tooLongHeaders < 20){ $tooLongHPos.= "$lineNum "};
						$warnings++;
					}
					if ($row =~ /^>\s.*/){
						$spaceNearName++;
						if ($spaceNearName < 20){$spaceNearNamePos.= "$lineNum "};
						$warnings++;
					}
			}else{
				#Here I chech if the sequence is too long
				if (length($row) > $suggestedLength){
					$tooLongSeqs++;
					if ($tooLongSeqs < 20){ $tooLongSPos.= "$lineNum "};
					$warnings++;
				}
				#Check if the sequence has IUPAC characters only
				if ($row !~ /^[$IUPACCode]*$/i){
					$noIUPACRows++;
					if ($noIUPACRows < 20){ $noIUPACRowsPos.= "$lineNum "};
					$warnings++;
				}	
			}			
	}
	my $nonzero = non_zero_seqs_length($fastaFilePath);
	if( $nonzero < 1){
		print "WARNING: The ".abs($nonzero)."-th sequence has length zero. You should not continue...\n";
		$warnings++;
	}
	#Here the user is alerted with a WARNING and row number
	if ( $tooLongHeaders > 0){
		print "WARNING: Found headers longer than $suggestedLength characters at lines: $tooLongHPos...\n";
	}
	if ( $tooLongSeqs > 0){
		print "WARNING: Found a sequence estending on the line more than $suggestedLength characters at lines: $tooLongSPos...\n";
	}
	if ( $spaceNearName > 0){
		print "WARNING: Found rows with a space after the symbol > at lines: $spaceNearNamePos...\n";
	}
	if ( $noIUPACRows > 0){
		print "WARNING: Found rows with non IUPAC characters at lines: $noIUPACRowsPos...\n";
		print "Allowed IUPAC characters can be found in the config_annocript.txt file.\n";
	}

	close(FASTA);
	
	if ($warnings == 0){
			print "The fasta file $fastaFilePath is good!\n";
	}
	return $warnings;
}


=head2 non_zero_seqs_length

 Title   : non_zero_seqs_length
 Usage   : non_zero_seqs_length( - configFilePath = path of the config file
                               );

 Function: checks if there are sequences of length zero in the fasta file.
					Means that the transcriptome is not well formed and you have
					two consecutive headers 
            
 Returns :  1 if the file has no zero length sequences
					-i if a sequence is found to have zero length and i is the number
					of sequence in the file

=cut
sub non_zero_seqs_length{  
	my $fasta = shift;
	
	my $seqin = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');
  my $seqNum = 1;                          
  my $found = 0;
  my $retVal = 0;
  while(my $seq = $seqin->next_seq) {
		if ($seq->length < 1){
			$found = 1;	
		}
		$seqNum++;
	}
	
	if ($found){
			$retVal = -($seqNum);
	}else{
			$retVal = 1;
	}
	
	return $retVal;
}

1;
