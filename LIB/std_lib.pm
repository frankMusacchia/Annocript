
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
    
    
#######################utilities.pm###########################################################
# std_lib.pm - A module that contains a series of utils subroutine                              #
################################################################################################

package LIB::std_lib;

BEGIN
{
use Exporter();

@ISA = qw(Exporter);
@EXPORT_OK = qw( print_array print_hash correct_type num_processors
				ram_memory save_hash deallocate_hash);
}

#GENERAL UTILITY MODULES
use strict;
use warnings;


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

=head2 print_hash

 Title   : print_hash
 Usage   : print_hash(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints the first two levels of an hash in input. 
 Returns : nothing

=cut
sub print_hashOLD{
  my $hash = shift;
  
	foreach my $key ( keys %{$hash}){
		#print $key.": ".$hash->{$key}."\n";			
		while (my ($intkey, $value) = each %{ $hash->{$key} } ) {
			print "$key->$intkey = ".$hash->{$key}->{$intkey}." \n";
		}
	}
}


=head2 print_hash

 Title   : print_hash
 Usage   : print_hash(  - array = the reference to the array to print (\@array)
                               );

 Function:  prints the first three level of an hash in input. 
 Returns : nothing

=cut
sub print_hash{
  my $hash = shift;
  
	foreach my $key ( keys %{$hash}){
		#print $key.": ".$hash->{$key}."\n";			
		while (my ($intkey, $value) = each %{ $hash->{$key} } ) {
			if (ref($value) eq 'HASH'){
					while (my ($intkey2, $value2) = each %{ $hash->{$key}->{$intkey}} ) {
						print "$key->$intkey-> $intkey2 = ".$hash->{$key}->{$intkey}->{$intkey2}." \n";
					}
			}else{
				print "$key->$intkey = ".$hash->{$key}->{$intkey}." \n";
			}
		}
	}
}

1;
