#!/usr/bin/perl

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

use strict;
use warnings;
use Cwd;
use File::Copy;


my $configUserFile = "CONFIGURATION/config_user.txt";

print "This script will prepare Annocript to work in your directory.\n";

# Are you sure you want to insert all here?
# If YES folders are created here and we can write a file for Annocript where are written the two paths of the folders.
# Finish!

my $sure= getInput("Are you sure you want to enjoy this software?(y or n)",'^[YynN]$');


if ( ($sure eq "y") or ($sure eq "Y")){ 
  
  #Asks to the user what is the directory where he wants to work
  my $workDir = '';    
  while ( !(-e $workDir) ){
    print "\nWrite an existing complete path where you want to play with Annocript (/home/username/ann_works): ";
    $workDir= <STDIN>;
    chomp $workDir; #to cancel the return 
  }
  
  #print "working dir: $workDir\n";
  #I take the folder path of Annocript
  my $annDir = getcwd;
  chmod 0755,"the_annocript.pl";
  
  #Writing a file with the paths to the current folder and the Annocript one
  open (FOLD, ">$workDir/folders.txt") or die "Cannot create file folders.txt. Exiting.. \n";
  print FOLD "$workDir $annDir";
  close(FOLD);

  #Moving the User configuration file in the working directory
  move($annDir."/$configUserFile",$workDir."/$configUserFile");
  print "Done! Now you can start Annocript from this folder!\n";
    
}else {die "Change directory and start again the install script.\n";}
  



=head2 getInput

 Title   : getInput
 Usage   : getInput( - sentence: a sentence that will printed in input to ask something to the user;
					 - regex: what the answer of the user have to respect
                               );

 Function: Takes in input a sentence and a regex. Asks to the user the sentence and controls its input with regex
 
 Returns : input given by the user

=cut
sub getInput{
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
