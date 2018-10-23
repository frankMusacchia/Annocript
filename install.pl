#!/usr/bin/perl

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

use strict;
use warnings;
use Cwd;
use File::Copy;
use Getopt::Long;

use LIB::mysql_db_management qw(check_DB_user_permissions );
#Using a library to manage files
use LIB::files_management qw( extract_name);
								
my $update = "";
								
my $configUserFile = "config_user.txt";

print "This script will prepare Annocript to work in your directory.\n";

# Are you sure you want to insert all here?
# If YES folders are created here and we can write a file for Annocript where are written the two paths of the folders.
# Finish!

my $sure= getInput("Are you sure you want to enjoy this software?(y or n)",'^[YynN]$');
my $d_programs_f = "DL_PROGRAMS";
my $install_dep_script = "install_all.sh";
my $programs_suff = " dna2pep* CPC* ";
my $sample_transcriptome = "trial_transcriptome.fasta";

parse_command_line_args();


if ( ($sure eq "y") or ($sure eq "Y")){ 
  
  #I take the folder path of Annocript
  my $annDir = getcwd;
  
  #Asks to the user what is the directory where he wants to work
  my $workDir = '';    
  while ( !(-e $workDir) ){
    print "\nWrite an existing complete path where you want to play with Annocript (/home/username/ann_works): ";
    $workDir= <STDIN>;
    chomp $workDir; #to cancel the return 
  }


	#################INSTALL DEPENDENCIES
	##This part of code allows the automated installation of all dependencies
	if ( $update > 0 ){ 
		print "The program will only be updated. Dependencies will not be installed.\n";
	}else{
		
		#Asks to the user what is the directory where he wants to install the programs
		my $binDir = '';    
		while ( !(-e $binDir) ){
			print "\nWrite an existing complete path where you want to install the programs (/home/username/bin/): ";
			$binDir= <STDIN>;
			chomp $binDir; #to cancel the return 
		}
		#Now go to the $d_programs_f folder
		chdir $d_programs_f;
				 
		#Install dependencies
		chmod 0755,$install_dep_script;
		system("./$install_dep_script $binDir/" ) == 1 or print "DONE..\n";

	}
	#################INSTALL DEPENDENCIES
	
  #print "working dir: $workDir\n";
	chdir $annDir;
  chmod 0755,"annocript.pl";
  chmod 0755,"USEFUL/faSomeRecords";
  
  #Writing a file with the paths to the current folder and the Annocript one
  open (FOLD, ">$workDir/folders.txt") or die "Cannot create file folders.txt. Exiting.. \n";
  print FOLD "$workDir $annDir";
  close(FOLD);

  #Copying the User configuration file in the working directory
  chmod 0755,$annDir."/$configUserFile";
  copy($annDir."/CONFIGURATION/$configUserFile",$workDir."/$configUserFile") or die "Unable to copy $annDir/CONFIGURATION/$configUserFile in $workDir/$configUserFile\n";
  
  #Copying the sample fasta transcriptome for the tutorial
  chmod 0755,$annDir."/$configUserFile";
  copy($annDir."/USEFUL/$sample_transcriptome",$workDir."/$sample_transcriptome") or die "Unable to copy $annDir/CONFIGURATION/$sample_transcriptome in $workDir/$sample_transcriptome\n";  
  print "Done! Now you can run Annocript from $workDir!\n";
    
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

=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.

	my $howToUse = "Use with: \nperl complete_install.pl \n\n".
	"\t-u|--update: Use this command if you just want to update VarGenius with the same dependencies you installed before\n";


	  
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "u|update" => \$update  );

	if ($update){$update=1;}else{$update=-1;}
  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

}
