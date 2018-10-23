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
use Data::Dumper;
use Time::HiRes qw( time );

#Annocript libraries
use LIB::db_creator qw( execute_db_creator);
use LIB::execute_programs qw( execute_programs run_blast_indexing);
use LIB::programs_out_parser qw(execute_programs_out_parser);
use LIB::output_and_stats qw(extractStatistics);


use LIB::files_management qw(print_file);


#The executor of Annocript is a perl script that can be used to run Annocript in background
#The Annocript configurator has previously checked the configuration file and built the configuration hashes
#then it created files with the stored hashes. 
#Therefore this script simply open those files and upload the hashes, then call the three Annocript modules
#to execute the db creation, programs and ouput writing.

my $configHash;
my $timesFile = $ARGV[2];
my $configHashPath = $ARGV[3];

#Words to match to be a non coding
my @trWords = ('uncharacterized','predicted protein','hypothetical','putative uncharacterized');

#open the timetable file
open(TIMES_FILE, ">$timesFile") or die "ERROR [$?]: Cannot open $timesFile! Check permissions.\n";

#Variables for computing execution times
my $globalStart = time;
my $partTime  = time;
my $partDuration = undef;


my $str1;
open(READ, "<$configHashPath") or die ("ERROR [$?]: Failure opening '$configHashPath \n"); 
{ local($/) = ""; $str1 = <READ>; close(READ); eval $str1; }#LOAD Hash FILE IN SCALAR

#Execution of the creation of the database
if ($configHash->{'doDbCreation'} eq "YES"){
  print "\n\n ####################### DB CREATION ########################\n";
  print "\n".scalar(localtime)." \nAnnocript is starting to create the database...\n";
  execute_db_creator($configHash);	
  print "\n ####################### ############ ########################\n";
  $partDuration = time-$partTime;
  print TIMES_FILE "dbcreation\t".$partDuration."\n";
  $partTime = time;
}

#Execution of programs 
if ($configHash->{'doExecutePrograms'} eq "YES"){
  print "\n\n ####################### PROGRAM EXECUTION ########################\n";
  print "\n".scalar(localtime)." \nAnnocript is starting the execution of programs...\n";
  execute_programs($configHash);
  print "\n ####################### ############ ########################\n";
  $partDuration = time-$partTime;
  print TIMES_FILE "exec_progs\t".$partDuration."\n";
  $partTime = time;
}

#Creation of the final table of output and gff3 files
if ($configHash->{'doBuildOutput'} eq "YES"){
  print "\n\n ####################### OUTPUT CREATION ########################\n";
  print "\n".scalar(localtime)." \nAnnocript is starting the generation of the results table...\n";
  execute_programs_out_parser($configHash);
  print "\n ####################### ############ ########################\n";
  $partDuration = time-$partTime;
  print TIMES_FILE "create_out\t".$partDuration."\n";
  $partTime = time;
}

#Extracting statistics and plots in a web page		
if ($configHash->{'extractStatistics'} eq "YES"){
  print "\n\n ####################### STATISTICS CREATION ########################\n";
  print "\n".scalar(localtime)." \nAnnocript is starting the generation of the statistics files and web page...\n";
	#statistics
	extractStatistics($configHash,@trWords);		
  print "\n ####################### ############ ########################\n";
  $partDuration = time-$partTime;
  print TIMES_FILE "create_stats\t".$partDuration."\n";
  $partTime = time;
}



#Calculate the time needed for all the processes
my $globalDuration = time - $globalStart;
print TIMES_FILE "total\t".$globalDuration."\n";
close (TIMES_FILE);

  
print "---->The Annocript is closed...............\n";

#Showing configuration files uses
print "\n\n\n\n CONFIGURATION FILES USED: \n\n";
print_file($ARGV[0]);
print "\n\n\n\n";
print_file($ARGV[1]);


