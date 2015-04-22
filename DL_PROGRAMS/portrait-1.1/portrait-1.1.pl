#!/usr/bin/perl
###########################################################################
#                                                                         #
#   This program is free software; you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 3 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program; if not, write to the                         #
#   Free Software Foundation, Inc.,                                       #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
#                                                                         #
###########################################################################
#                                                                         #
#   PORTRAIT Prediction of transcriptomic ncRNA by ab initio methods      #                                                        #                                                                         #
#                                                                         #
#   Description:                                                          #
#      an algorithm for identifying non-coding RNAs on a transcriptome    #
#      context.                                                           #
#                                                                         #
#   Developed by:                                                         #
#                                                                         #
###########################################################################
use strict;
use Getopt::Long;
$| = 1;

#RTARRIAL Jun/2007~Feb/2008

##########################################################################
#  PORTRAIT -- Prediction of Transcriptome ncRNA by Ab Initio Methods
##########################################################################
# Utility: given a input file containing transcritome (nucleotide) sequences,
#   PORTRAIT returns probability of each sequence being a ncRNA.
# Requirements: Currently we support only UNIX-based OS.
#   Programs required: CAST 1.0, ANGLE and LIBSVM-2.84. All these third party
#   programs should be included in the package. LIBSVM is the only one
#   that requires installation.
#   PORTRAIT also requires two other PERL scripts and two SVM models,
#   all should also be included in the package.
##############
# Input format: Currently we support only FASTA-formatted, nucleotide sequences.
#   Your sequences may only contain these characters: [ATGCUN]. If a sequence contain
#   other characters this sequence will be ignored.
#   Also notice that sequences with more than 20% of their composition composed of 'N'
#   will be ignored. This is because too many 'N's may degrade statistical analysis performance.
#   Sequences with less than 80 letters will also be ignored. If length is less than 80nt,
#   ORF prediction is innacurate, and your sequence is very likely non coding.
#   If any sequences are ignored, a log file will be created reporting the errors.

# Please report any bugs or feedback to rtarrial[at]gmail.com

# Don´t modify anything in this script, unless you know what you are doing.

### INITIAL CONFIG - NO NEED TO CHANGE ###
my $home_dir = '/home/rtarrial/prog'; # Address of this script
my $libsvm_dir = '/home/rtarrial/libsvm-2.84'; # Third-party
my $cast_dir = '/home/rtarrial/cast-linux'; # Third-party
my $angle_dir = '/home/rtarrial/angle'; # Third-party
my $wrapper_script_ptn = 'SVM_wrapper_2_dual.pl'; # NO NEED TO CHANGE! calculates attributes and wraps them for transcripts with orfs
my $wrapper_script_nt = 'SVM_wrapper_2_nt.pl';    # NO NEED TO CHANGE! same as above, but for transcripts without orfs, so only nucleotide attributes are determined
my $ptn_model = 'withProteins.model';        #No need to change!
my $nt_model = 'noOrf.model';  #No need to change!
my $config_file = 'PORTRAIT.config'; # Internal Config file.
my $angle_bit = '32'; #this may be changed in the config process; don´t worry about it
######## END USER CONFIG ############
## NO USER SERVICEABLE PARTS AHEAD ##

#### CONFIGURATION PROCESS ####
if (-e $config_file) {
 open(CONFIG, "$config_file") or die "\nUnable to open $config_file: $!\n";
 my $right_config = 0;
 my $count_progs = 0;
 READCONFIG:while (<CONFIG>) {
  chomp $_;
  if ($_ =~ /^#PORTRAIT CONFIG FILE/) {$right_config = 1; next READCONFIG;}
  elsif (($_ =~ /^home::(.*)/) && ($right_config)) {$home_dir = $1; $count_progs++;}
  elsif (($_ =~ /^libsvm::(.*)/) && ($right_config)) {$libsvm_dir = $1; $count_progs++;}
  elsif (($_ =~ /^cast::(.*)/) && ($right_config)) {$cast_dir = $1; $count_progs++;}
  elsif (($_ =~ /^angle::(.*)/) && ($right_config)) {$angle_dir = $1; $count_progs++;}
  elsif (($_ =~ /^angle_bit::(.*)/) && ($right_config)) {$angle_bit = $1; $count_progs++;}
 }
 if ($count_progs != '5') {unlink ('$config_file'); die "\nInvalid config file. Please run script again.\n";}
}
else {
 print "\nWelcome! Now starting configuration process.";
 print "\nDirectory where PORTRAIT is located (e.g. /home/user/program/): ";
 $home_dir = <STDIN>; chomp $home_dir;
 if ($home_dir =~ /(.*)\/$/) {$home_dir = $1;} #get rid of last slash, if any
 unless (-e "$home_dir/$wrapper_script_ptn") {die "\nERROR: $wrapper_script_ptn not found in specified dir!\n";}
 unless (-e "$home_dir/$wrapper_script_nt") {die "\nERROR: $wrapper_script_ptn not found in specified dir!\n";}
 unless (-e "$home_dir/$0") {die "\nERROR: $0 not found in specified dir!\n";}
 unless (-e "$home_dir/$ptn_model") {die "\nERROR: $ptn_model not found in specified dir!\n";}
 unless (-e "$home_dir/$nt_model") {die "\nERROR: $nt_model not found in specified dir!\n";}
 print "\nDirectory where LIBSVM is installed (e.g. /home/user/libsvm-2.84/): ";
 $libsvm_dir = <STDIN>; chomp $libsvm_dir;
 if ($libsvm_dir =~ /(.*)\/$/) {$libsvm_dir = $1;}
 unless (-e "$libsvm_dir/svm-predict") {die "\nERROR: svm-predict not found in specified dir!\n";}
 testLibsvm($libsvm_dir);
 print "\nDirectory where CAST is installed (e.g. /home/user/cast-linux/): ";
 $cast_dir = <STDIN>; chomp $cast_dir;
 if ($cast_dir =~ /(.*)\/$/) {$cast_dir = $1;}
 unless (-e "$cast_dir/cast"){die "\nERROR: cast not found in specified dir!\n";}
 testCast($cast_dir);
 print "\nDirectory where ANGLE is installed (e.g. /home/user/angle/): ";
 $angle_dir = <STDIN>; chomp $angle_dir;
 if ($angle_dir =~ /(.*)\/$/) {$angle_dir = $1;}
 unless (-e "$angle_dir/ANGLE-linux32DP"){die "\nERROR: ANGLE 32 not found in specified dir!\n";}
 unless (-e "$angle_dir/ANGLE-linux64DP"){die "\nERROR: ANGLE 64 not found in specified dir!\n";}
 $angle_bit = testAngle($angle_dir);
 print "Running ANGLE ${angle_bit}-bit version.\n\n";
 open(CONFIG, ">$config_file") or die "\nConfig file: $!\n";
 print CONFIG "#PORTRAIT CONFIG FILE\nhome::$home_dir\nlibsvm::$libsvm_dir\ncast::$cast_dir\nangle::$angle_dir\nangle_bit::$angle_bit";
 die "Config process successful! Run again to start working.\n";
}
######### END CONFIG PROCESS #########


sub usage {
    my $prog_name = $0;
    print <<END;
\n./$prog_name -i <transcripts file> -s <libsvm> -c <cast> -a <angle> [options]
	-i <transcripts file> : input file with multiple FASTA-formatted NUCLEOTIDE sequences.
        -s <LIBSVM path> : OPTIONAL provide path where LIBSVM software is located.
        -c <CAST path> : OPTIONAL provide path where CAST software is located.
        -a <ANGLE path> : OPTIONAL provide path where ANGLE software is located.
        -b [32|64] : OPTIONAL whether to use ANGLE 32- or 64-bit version.
        -t [0|1] : OPTIONAL toggle deletion of temporary files.
		1 : temp files are kept. (Used to speed up subsequent runs on same test set and debugging)
		0 : temp files are deleted. [DEFAULT]
        -q [0|1] : OPTIONAL toggle Quiet mode.
                0 : verbose mode. Progress is reported on STDOUT. [DEFAULT]
                1 : quiet mode. Only error messages and LIBSVM output are reported.
        -h : This help.\n
END
    exit;
}

my $input_transcript_original_file = '&';
my $transcripts_with_orfs;
my $transcripts_without_orfs;
my $attrib_mode;
my $keep_temp_files = 0;
my $ask_for_help = 0;
my $quiet = 0;

usage() if ($#ARGV<0);

Getopt::Long::Configure ("bundling", "no_ignorecase");
GetOptions (
    "t" => \$keep_temp_files,
    "h" => \$ask_for_help,
    "q" => \$quiet,
    "b=s" => \$angle_bit,
    "i=s" => \$input_transcript_original_file,
    "s=s" => \$libsvm_dir,
    "c=s" => \$cast_dir,
    "a=s" => \$angle_dir,
);

usage() if $ask_for_help;

#getting rid of last slash
if ($libsvm_dir =~ /(.*)\/$/) {$libsvm_dir = $1;}
if ($cast_dir =~ /(.*)\/$/) {$cast_dir = $1;}
if ($angle_dir =~ /(.*)\/$/) {$angle_dir = $1;}

if ($input_transcript_original_file eq '&') {die "\nField -i is mandatory!\n";}
unless (($angle_bit == '32') or ($angle_bit == '64')) {die "\nField -b must be \"32\" or \"64\" -- invalid designation: $angle_bit\n";}

#This check was already done when config file was created, but to make sure all programs are there, this redundant check is made

unless ((-e "${libsvm_dir}\/svm-predict") && ("${libsvm_dir}\/svm-train")) {die "\nLIBSVM not found! Provided: $libsvm_dir. Example: \"\/usr\/bin\/libsvm-2.84\"\n";}
testLibsvm($libsvm_dir);
unless (-e "${cast_dir}\/cast") {die "\nCAST not found! Provided: $cast_dir. Example: \"\/usr\/bin\/cast-linux\"\n";}
testCast($cast_dir);
unless (-e "$angle_dir") {die "\nANGLE not found! Provided: $angle_dir. Example: \"\/usr\/bin\/angle\"\n";}
testAngleSecond($angle_dir,$angle_bit);
existe_arquivo ($input_transcript_original_file);
my $log_file = "$input_transcript_original_file.log";
if (-e $log_file) {unless($quiet){print "Log file found: overwriting!\n";}}
my $transcript_original_file = "${input_transcript_original_file}_results";
nao_existe_arquivo ($transcript_original_file);
open(COPY_INPUT, "$input_transcript_original_file") or die "\nError opening input $transcript_original_file: $!\n";
open(COPIED_INPUT, ">$transcript_original_file") or die "\nError creating copy $transcript_original_file: $! $!\n";
while (<COPY_INPUT>) {
 print COPIED_INPUT $_;
}
close COPY_INPUT;
close COPIED_INPUT;
#There, now we dont need input file from user and may modify input file as we like


##### DOS2UNIX #####
# This *SHOULD* work. Anyway, we recommend you to use dos2unix in your sequence before submission
system("perl -pi -e 's/\r\n/\n/g;' $transcript_original_file");
#####
fix_input($transcript_original_file, $log_file);

$transcripts_with_orfs = "${transcript_original_file}.withorfs";
#nao_existe_arquivo($transcripts_with_orfs);
$transcripts_without_orfs = "${transcript_original_file}.withoutorfs";
my $orf_file = "${transcript_original_file}.orf";
#nao_existe_arquivo($transcripts_without_orfs);
nao_existe_arquivo("${transcript_original_file}.SVMpredictions");
nao_existe_arquivo("${transcript_original_file}_noOrf.SVMpredictions");
nao_existe_arquivo("${transcript_original_file}_ncRNA.scores");
nao_existe_arquivo("${transcript_original_file}_noOrf_ncRNA.scores");
nao_existe_arquivo("${transcript_original_file}_noOrf_ncRNA.fas");
nao_existe_arquivo("${transcript_original_file}_ncRNA.fas");

my @current_dir_contents = <*>;
my $ptn_wrapper_ok = 0;
my $nt_wrapper_ok = 0;
my $ptn_model_ok = 0;
my $nt_model_ok = 0;
foreach my $file (@current_dir_contents) {
  if ($file eq "$wrapper_script_ptn") {$ptn_wrapper_ok = '1';}
  elsif ($file eq "$wrapper_script_nt") {$nt_wrapper_ok = '1';} 
  elsif ($file eq "$ptn_model") {$ptn_model_ok = '1';}
  elsif ($file eq "$nt_model") {$nt_model_ok = '1';}
}
unless ($ptn_wrapper_ok) {die "\nFATAL :: Could not locate script $wrapper_script_ptn in same dir as $0!\n";}
unless ($nt_wrapper_ok) {die "\nFATAL :: Could not locate script $wrapper_script_nt in same dir as $0!\n";}
unless ($nt_model_ok) {die "\nFATAL :: Could not locate script $ptn_model in same dir as $0!\n";}
unless ($nt_model_ok) {die "\nFATAL :: Could not locate script $nt_model in same dir as $0!\n";}

my $dir_with_files = $home_dir;
unless (-e "$dir_with_files/$wrapper_script_ptn") {die "\nFATAL :: Perl script \"$wrapper_script_ptn\" must be in same Directory as $0!\n";}
unless (-e "$dir_with_files/$wrapper_script_nt") {die "\nFATAL :: Perl script \"$wrapper_script_nt\" must be in same Directory as $0!\n";}
unless (-e "$dir_with_files/$ptn_model") {die "\nFATAL :: SVM model file \"$ptn_model\" must be in same Directory as $0!\n";}
unless (-e "$dir_with_files/$nt_model") {die "\nFATAL :: SVM model file \"$nt_model\" must be in same Directory as $0!\n";}

my $time_now = gettime();
unless ($quiet) { print "\nStart time: $time_now\n----------------------";}
unless ($quiet) {print "\nTranslating transcripts...";}
if (-e $orf_file) { unless ($quiet) {print "file found, skip.";}}
else {
 angle_translation ($transcript_original_file, $orf_file, $transcripts_with_orfs, $transcripts_without_orfs, $angle_dir, $angle_bit);
 unless ($quiet) {print "OK.";}
}

#Checking transcript status. Will skip protein processing if
#no orfs were found in the input

 my $count_withorfs = 0;
 my $count_withoutorfs = 0;
 open(WITHORFS, $transcripts_with_orfs) or die "\n$!\n";
 while (<WITHORFS>) {
  if ($_ =~ /^>/) {$count_withorfs++;}
 }
 close WITHORFS;
 open(WITHOUTORFS, $transcripts_without_orfs) or die "\n$!\n";
 while (<WITHOUTORFS>) {
  if ($_ =~ /^>/) {$count_withoutorfs++;}
 }
 close WITHOUTORFS;
 if (($count_withorfs == 0) && ($count_withoutorfs == 0)) { #the "withoutorfs" file is generated by angle translation, ie, no withorfs = no withoutorfs
 # $transcripts_without_orfs = $transcript_original_file;
   open(WITHOUTORFS, ">$transcripts_without_orfs") or die "\n$!\n";
   open(ORIGINAL, "$transcript_original_file") or die "\n$!\n";
   select WITHOUTORFS; $|=1; select ORIGINAL; $|=1; select STDOUT;
   while(<ORIGINAL>) {
    print WITHOUTORFS $_;
   }
  close WITHOUTORFS;
  close ORIGINAL;
  open(WITHOUTORFS, "$transcripts_without_orfs") or die "\n$!\n";
   while (<WITHOUTORFS>) {
    if ($_ =~ /^>/) {$count_withoutorfs++;}
   }
   close WITHOUTORFS;
 }
 unless ($quiet) {print "\nFound $count_withorfs transcripts with ORFs and $count_withoutorfs without ORFs.";}
 my $all_with_orfs = 0;
 my $all_without_orfs = 0;
 if ($count_withoutorfs == 0) {$all_with_orfs = 1;}
 elsif ($count_withorfs == 0) {$all_without_orfs = 1;}



## BEGIN SANITY CHECK ##
my $count_transcripts = 0;
my $count_proteins = 0;
unless ($all_with_orfs) {$count_transcripts = checkFastaNT($transcripts_with_orfs);}
unless ($all_without_orfs) {$count_proteins = checkFastaPT($orf_file);}
unless (($all_with_orfs) or ($all_without_orfs)) {
 if ($count_transcripts != $count_proteins) {die "\nFATAL :: Amount of ORFs ($count_proteins) and corresponding transcripts ($count_transcripts) differ. Please check again your files.\n";}
}
## END SANITY CHECK ##

 my $iep_raw_output;
 my $cast_raw_output;
 my $iep_parsed_output;
 my $cast_parsed_output;

unless ($all_without_orfs) {
 $iep_raw_output = "${orf_file}.iep_raw";
 $cast_raw_output = "${orf_file}.cast_raw";
 $iep_parsed_output = "${orf_file}.iep_parsed";
 $cast_parsed_output = "${orf_file}.cast_parsed";
}

#
unless ($quiet) {print "\nGathering and wrapping attributes...";} #An error at this point may indicate dos formatting problem. Solve this using dos2unix utility
unless ($all_without_orfs) {
 unless (-e $cast_parsed_output) {
  execute_cast($orf_file, $cast_raw_output, $cast_dir);
  parse_cast($cast_raw_output, $cast_parsed_output);
 }
}
#############
# Some atributes collected, we now gather other, join and format them calling another script
unless ($all_without_orfs) {
 system("perl $dir_with_files/$wrapper_script_ptn $transcripts_with_orfs $orf_file 0") == 0 or die "\nUnable to call or error calling script $wrapper_script_ptn\n";
 unless ($quiet) {print "\nProtein-dependent OK.";}
}
unless ($all_with_orfs) {
	my $command = "perl $dir_with_files/$wrapper_script_nt $transcripts_without_orfs 0";#FRANKCODE
	
 system($command) == 0 or die "\nUnable to execute command $command\n";#FRANKCODE
 unless ($quiet) {print "\nNucleotide-only OK.";}
}
unless ($keep_temp_files) {
 unlink(
# "$orf_file", # Might be interesting to see the orf file
 "${transcripts_with_orfs}_N.rnasize",
 "${transcripts_with_orfs}_N.nt1",
 "${transcripts_with_orfs}_N.nt2",
 "${transcripts_with_orfs}_N.nt3",
 "${orf_file}_P.orfsize",
 "${orf_file}_P.aacomp",
 "${orf_file}_P.soap",
 "${orf_file}_P.iep",
 "$cast_parsed_output",
 "$cast_raw_output",
 "${transcripts_without_orfs}_N.rnasize",
 "${transcripts_without_orfs}_N.nt1",
 "${transcripts_without_orfs}_N.nt2",
 "${transcripts_without_orfs}_N.nt3",
 );
}
#############
# Run LIBSVM and gather predictions
my $trash = 'undef';
unless ($quiet) {print "\nRunning LIBSVM...\n";}
unless ($all_without_orfs) {
 unless ($quiet) {system("$libsvm_dir/svm-predict -b 1 ${transcripts_with_orfs}.SVM ${dir_with_files}/$ptn_model ${transcript_original_file}.SVMpredictions") == 0 or die "\nUnable to run LIBSVM step 1: $!\n";}
 else {$trash = `$libsvm_dir/svm-predict -b 1 ${transcripts_with_orfs}.SVM ${dir_with_files}/$ptn_model ${transcript_original_file}.SVMpredictions 2>&1`; if ($trash eq 'undef') {die "\nUnable to run LIBSVM step 1: $!\n";}}
}
unless ($all_with_orfs) {
 unless ($quiet) {system("$libsvm_dir/svm-predict -b 1 ${transcripts_without_orfs}.SVM ${dir_with_files}/$nt_model ${transcript_original_file}_noOrf.SVMpredictions") == 0 or die "\nUnable to run LIBSVM step 2: $!\n";}
 else {$trash = `$libsvm_dir/svm-predict -b 1 ${transcripts_without_orfs}.SVM ${dir_with_files}/$nt_model ${transcript_original_file}_noOrf.SVMpredictions 2>&1`; if ($trash eq 'undef') {die "\nUnable to run LIBSVM step 2: $!\n";}}
}
unless ($quiet) {print "OK.\n";}
$trash = '';
############
# Grab headers of transcripts predicted as ncRNA
unless ($all_without_orfs) {
 generate_ncrna_headers("${transcript_original_file}.SVMpredictions", $transcripts_with_orfs, "${transcript_original_file}_ncRNAs.scores", "${transcript_original_file}_ncRNA.fas");
}
unless ($all_with_orfs) {
 generate_ncrna_headers("${transcript_original_file}_noOrf.SVMpredictions", $transcripts_without_orfs, "${transcript_original_file}_noOrf_ncRNA.scores", "${transcript_original_file}_noOrf_ncRNA.fas");
}

if ($all_without_orfs) {
 open(CREATE_COD, ">${transcript_original_file}_ncRNA.fas") or die "\n$!\n";
 close CREATE_COD;
 open(SCORES_COD, ">${transcript_original_file}_ncRNAs.scores") or die "\n$!\n";
 close SCORES_COD;
}
elsif ($all_with_orfs) {
 open(CREATE_NC, ">${transcript_original_file}_noOrf_ncRNA.fas") or die "\n$!\n";
 close CREATE_NC;
 open(SCORES_NC, ">${transcript_original_file}_noOrf_ncRNA.scores") or die "\n$!\n";
 close SCORES_NC;
}

 #We join ORF-independent and ORF-dependent predicted ncRNAs
 open(RECEPTOR, ">>${transcript_original_file}_ncRNA.fas") or die "\n$!\n";
 open(DONOR, "${transcript_original_file}_noOrf_ncRNA.fas") or die "\n$!\n";
 while (<DONOR>) {
  print RECEPTOR $_;
 }
 close DONOR;
 close RECEPTOR;
 open (DONOR1, "${transcript_original_file}_ncRNAs.scores") or die "\n$!\n";
 open (DONOR2, "${transcript_original_file}_noOrf_ncRNA.scores") or die "\n$!\n";
 open (RECEPTOR, ">${transcript_original_file}_ncRNA.scores") or die "\n$!\n";
 print RECEPTOR "[SEQUENCE] ::: [Coding Probability] [Non-Coding Probability]\n";
 while (<DONOR1>) {
  print RECEPTOR $_;
 }
 while (<DONOR2>) {
  print RECEPTOR $_;
 }
 close DONOR;
 close RECEPTOR;

 my $all_scores = "${transcript_original_file}_all.scores";

 if ($all_with_orfs) {
  open(ALL, ">$all_scores") or die "\n$all_scores : $!\n";
  open(PTN, "${transcript_original_file}_ncRNAs.scores") or die "\n${transcript_original_file}_ncRNAs.scores : $!\n";
  while (<PTN>) {
   print ALL $_;
  }
 close ALL;
 close PTN;
 convert_scores($all_scores);
 }
 elsif ($all_without_orfs) {
  open(ALL, ">$all_scores") or die "\n$all_scores : $!\n";
  open(NT, "${transcript_original_file}_noOrf_ncRNA.scores") or die "\n${transcript_original_file}_noOrf_ncRNA.scores : $!\n";
  while (<NT>) {
   print ALL $_;
  }
 close ALL;
 close NT;
 convert_scores($all_scores);
 }
 else {
  generate_general_headers ("${transcript_original_file}.SVMpredictions", $transcripts_with_orfs, "${transcript_original_file}_ptn_all.scores");
  generate_general_headers ("${transcript_original_file}_noOrf.SVMpredictions", $transcripts_without_orfs, "${transcript_original_file}_nt_all.scores");
 }

unless (($all_with_orfs) or ($all_without_orfs)) {
 open(ALL, ">$all_scores") or die "\n$all_scores : $!\n";
 open(PTN, "${transcript_original_file}_ptn_all.scores") or die "\n${transcript_original_file}_ptn_all.scores : $!\n";
 open(NT, "${transcript_original_file}_nt_all.scores") or die "\n${transcript_original_file}_nt_all.scores : $!\n";

 while (<PTN>) {print ALL $_;}
 while (<NT>) {print ALL $_;}

 close ALL;
 close PTN;
 close NT;
}

 unless ($keep_temp_files) {
 unlink(
 "${transcripts_with_orfs}.SVM",
 "${transcripts_without_orfs}.SVM",
 "$transcript_original_file",
 "${transcript_original_file}_noOrf_ncRNA.fas",
 "${transcripts_with_orfs}",
 "${transcripts_without_orfs}",
 "${transcript_original_file}.SVMpredictions",
 "${transcript_original_file}_noOrf.SVMpredictions",
 "${transcript_original_file}_noOrf_ncRNA.scores",
 "${transcript_original_file}_noOrf.SVMpredictions",
 "${transcript_original_file}.SVMpredictions",
 "${transcript_original_file}_ncRNA.scores",
 "${transcript_original_file}_noOrf_ncRNAs.scores",
 "${transcript_original_file}_ncRNAs.scores",
 "${transcript_original_file}_noOrf_ncRNA.fas",
 "${transcript_original_file}_ptn_all.scores",
 "${transcript_original_file}_nt_all.scores"
 );
 } #end if keep temp files

unless ($quiet) {print "\nORFs file : $orf_file\nFASTA file  : ${transcript_original_file}_ncRNA.fas\nScores file (all) : ${transcript_original_file}_all.scores\n";}
$time_now = gettime();
unless ($quiet) {print "----------------------\nEnd time: $time_now\n\n";}

##################################
#      SUBROTINAS             #
##################################
sub existe_arquivo {
        my ($arquivo) = @_;
        unless (-e $arquivo) { die "\nFATAL :: File $arquivo does not exist!\n"; }
}

##################################

sub nao_existe_arquivo {
        my ($arquivo) = @_;
        if (-e $arquivo) { die "\nFATAL :: Output file $arquivo already exists!\n";} else { ; }
}

##################################
sub existe_dir {
  my ($pasta) = @_;
  unless (-d $pasta) {die "\nFATAL :: Specified dir $pasta does not exist!\n";}
}

##################################

sub abre_arquivo {
        my ($arquivo) = @_;
        my @conteudo = ();
        open (ARQ, $arquivo) or die "\nFATAL :: (sub abre_arquivo)::Could not open $arquivo.\n";
        select ARQ; $| = 1; select STDOUT;
        @conteudo = <ARQ>;
        close ARQ;
        return @conteudo;
}
##############################

sub separate_transcripts_without_orfs {
my ($in_fasta, $orfpred, $output, $output_withoutorfs) = @_;

my @pred_heads = ();
my @conteudo_pred = abre_arquivo($orfpred);
foreach my $head_pred (@conteudo_pred) {
 if ($head_pred =~ /^>/) {
  my @get_id = split (/\t/, $head_pred); #pegar somente a id, dispensar outras infos do arquivo com orfs
  push (@pred_heads, $get_id[0]);
 }
}
unless ($quiet) {print "\nComparing ORFs and original FASTA...";}
@conteudo_pred = ();
open(INFASTA, $in_fasta) or die "\nCould not open $in_fasta: $!\n";
open(OUTPUT, ">$output") or die "\nCould not save to $output: $!\n";
open(OUTPUT_WITHOUTORFS, ">$output_withoutorfs") or die "\nCould not save to $output_withoutorfs: $!\n";
select INFASTA; $| = 1; select OUTPUT; $| = 1; select OUTPUT_WITHOUTORFS; $| = 1; select STDOUT;

my $h = '0';
my $sequencia;
my $header;
my $notfirsttime = '0';
my $linha;

while ($linha = <INFASTA>) {
  if ($h > $#pred_heads) { close INFASTA;}
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
    if ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco

  if (($linha =~ /^>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
   chomp $header;
   if ($header eq $pred_heads[$h]) {
     print OUTPUT "$header\n";
     print OUTPUT "$sequencia\n"; #por motivo desconhecido a sequencia nao eh imprimida se estiver junto do header
     $h++;
   }
  else {print OUTPUT_WITHOUTORFS "$header\n"; print OUTPUT_WITHOUTORFS "$sequencia\n";}
   $sequencia = '';
   $header = $linha;
  }
  elsif(($linha =~ /^>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $header = $linha;
    $notfirsttime = '1';
  }
  elsif($linha !~ /^>/) { #poderia ser else
   $sequencia .= $linha;
   chomp $sequencia;
   $sequencia =~ s/\s//g;
   $sequencia =~ s/\n//g;
  } #Processando ultima seq

} #end while - arquivo fasta original
   chomp $header;
   if ($header eq $pred_heads[$h]) { #In "use warnings" mode, we may be warned that $pred_heads[$h] does not exist. Just ignore this message.
      print OUTPUT "$header\n";
      print OUTPUT "$sequencia\n";
   }
close INFASTA;
unless ($quiet) {print "OK.";}
} #end sub separate_transcripts_without_orfs
####################################
sub execute_cast {
 my ($input, $output, $cast_dir) = @_;
 if (-e $output) {return;}
 else {
  system("$cast_dir/cast $input > $output");
 }
} #end sub execute cast
####################################
sub parse_cast {
        my ($arquivo, $output) = @_;
        if (-e $output) {return;}
        my $notfirsttime = '0';
        my $sequencia;
        my $X = '0';

 open(OUTPUT, ">$output") or die "\n$output -- $!\n";
 select OUTPUT; $|=1; select STDOUT;
 open(ARQUIVO, $arquivo) or die "\n$arquivo -- $!\n";
        while (my $linha = <ARQUIVO>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  if ($linha =~ /^\s*$/) {next;} #ignorando campos em branco

  if (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
    $X = ($sequencia =~ tr/Xx//);
    $X = sprintf("%.17f", $X/length($sequencia));
    print OUTPUT "$X\n";
    $sequencia = '';
  }
  elsif(($linha =~ /^\>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $notfirsttime = '1';
  }
  elsif($linha !~ /^\>/) {
   $sequencia .= "${linha}";
   chomp $sequencia;
   $sequencia =~ s/\s//;
   $sequencia =~ s/\n//;
  }

} #final foreach; processando ultima seq

    $X = ($sequencia =~ tr/Xx//);
    $X = sprintf("%.17f", $X/length($sequencia));
    print OUTPUT "$X\n";
close ARQUIVO;
close OUTPUT;
return;
} #end sub parse cast
####################################

sub generate_ncrna_headers {
 my ($preds_file, $fasta_file, $output_file, $output_fasta) = @_;
 my @preds = ();
 my @confianca_preds = ();
 my $count = 0;
open(PREDS, $preds_file) or die "\n$preds_file -- $!\n";
 LOOPRED:while (my $line_preds = <PREDS>) {
  chomp $line_preds;
  if ($line_preds eq 'labels 1 -1') {next LOOPRED;} #nao iniciar counter!
  $count++;
  if ($line_preds =~ /-1 (.*)/) {
   push (@preds, $count); push (@confianca_preds, $1);
  }
 }
 if ($#preds != $#confianca_preds) {die "\nError: amount of preds and confidence differ?\n";}
 close PREDS;
#print @preds;
#Lista de posicoes esta correta

my $count_fasta = 0;

existe_arquivo($fasta_file);
nao_existe_arquivo($output_file);

my $notfirsttime = '0';
my $sequencia;
my $header;
my $confpred;

open(FASTA, $fasta_file) or die "\n$!";
open(OUT_FASTA, ">$output_fasta") or die "\n$output_fasta -- $!\n";
open(OUT_SCORES, ">$output_file") or die "\n$output_file -- $!\n";
select FASTA; $| = 1; select OUT_FASTA; $| = 1; select OUT_SCORES; $| = 1; select STDOUT;
my $advance_pred = 0;

ARQLOOP:while (my $linha = <FASTA>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  elsif ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco
##################---HEAD---######################
  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
   $count_fasta++;
   if ($advance_pred > $#preds) {last ARQLOOP;}
   if ($count_fasta == $preds[$advance_pred]) {
      chomp $header;
      $confpred = shift @confianca_preds;
      print OUT_SCORES "$header ";
      print OUT_SCORES "::: $confpred\n";
      print OUT_FASTA "$header\n";
      print OUT_FASTA "$sequencia\n";
      $advance_pred++;
   }
  $sequencia = '';
  $header = $linha;
  }
##################---HEAD---######################
  elsif(($linha =~ /^\>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $header = $linha;
    $notfirsttime = '1';
  }
  elsif($linha !~ /^\>/) { #poderia ser else
   $sequencia .= $linha;
   chomp $sequencia;
   $sequencia =~ s/\s//g;
   $sequencia =~ s/\n//g;
  }

} #final while; processando ultima seq
##################---HEAD---######################
   $count_fasta++;
#   if ($count_fasta <= $#preds) {
    if ($count_fasta == $preds[$advance_pred]) {
       chomp $header;
       $confpred = shift @confianca_preds;
       print OUT_SCORES "$header ";
       print OUT_SCORES "::: $confpred\n";
       print OUT_FASTA "$header\n";
       print OUT_FASTA "$sequencia\n";
       $advance_pred++;
    }
#   }
close FASTA;
close PREDS;
close OUT_SCORES;
close OUT_FASTA;
return;
} #end sub generate ncrna headers

#############################
sub generate_general_headers {
 my ($preds_file, $fasta_file, $output_file) = @_;
 my @preds = ();
 my @confianca_preds = ();
 my $count = 0;

existe_arquivo($fasta_file);
nao_existe_arquivo($output_file);

open(PREDS, $preds_file) or die "\n$preds_file -- $!\n";
 LOOPRED:while (my $line_preds = <PREDS>) {
  chomp $line_preds;
  if ($line_preds eq 'labels 1 -1') {next LOOPRED;} #nao iniciar counter!
  else {push (@confianca_preds, $line_preds);}
  }

 close PREDS;
#print @preds;
#Lista de posicoes esta correta


my $notfirsttime = '0';
my $sequencia;
my $header;
my $confpred;

open(FASTA, $fasta_file) or die "\n$!";
open(OUT_SCORES, ">$output_file") or die "\n$output_file -- $!\n";
select FASTA; $| = 1; select OUT_SCORES; $| = 1; select STDOUT;
my $class;
my @split_conf;
my $count = 0;

ARQLOOP:while (my $linha = <FASTA>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  elsif ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco
##################---HEAD---######################
  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
      chomp $header;
      $confpred = $confianca_preds[$count];
      $count++;
      @split_conf = ();
      @split_conf = split(/ /,$confpred);
      if ($split_conf[0] eq '1') {$class = 1;}
      elsif ($split_conf[0] eq '-1') {$class = 0;}
      print OUT_SCORES "$header ";
      print OUT_SCORES ":$class: $split_conf[1] $split_conf[2]\n";
  $sequencia = '';
  $header = $linha;
  }
##################---HEAD---######################
  elsif(($linha =~ /^\>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $header = $linha;
    $notfirsttime = '1';
  }
  elsif($linha !~ /^\>/) { #poderia ser else
   $sequencia .= $linha;
   chomp $sequencia;
   $sequencia =~ s/\s//g;
   $sequencia =~ s/\n//g;
  }

} #final while; processando ultima seq
##################---HEAD---######################
      chomp $header;
      $confpred = $confianca_preds[$count];
      @split_conf = ();
      @split_conf = split(/ /,$confpred);
      if ($split_conf[0] eq '1') {$class = 1;}
      elsif ($split_conf[0] eq '-1') {$class = 0;}
      print OUT_SCORES "$header ";
      print OUT_SCORES ":$class: $split_conf[1] $split_conf[2]\n";

#   }
close FASTA;
close PREDS;
close OUT_SCORES;
close OUT_FASTA;
return;
} #end sub generate general headers

############################
sub convert_scores {
 my ($file) = @_;
 open (IN, $file) or die "\nConverting score: $file: $!\n";
 my $temp_file = "${file}_convert.temp";
 open(TEMP, ">$temp_file") or die "\nConverting score: $temp_file: $!\n";
 my $class;
 while (<IN>) {
  chomp $_;
  if ($_ =~ /(^>.*) ::: (\S*) (\S*)/) {
   my $field1 = $1; my $field2 = $2; my $field3 = $3;
      if ($field3 > 0.5) {$class = 0;}
      elsif ($field2 > 0.5) {$class = 1;}
      print TEMP "$field1 ";
      print TEMP ":$class: $field2 $field3\n";
  }
 }
 unlink $file;
 rename($temp_file, $file);
} #end sub convert_scores

############################
sub gettime {
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900; #Y2K-compliant
return ("$hour\:$min\:$sec - $mday\/$mon\/$year");
} #end subrout gettime
#############################
#This frontend receives you multiFASTA, split it, then runs and parses ANGLE output, and formats it as protein.
#Special version with DP score output and unlimited sequence length acceptance.
#Software used with permission from Dr. Kana Shimizu. Reference:
#SHIMIZU, K., ADACHI, J., MURAOKA, Y. ANGLE: A SEQUENCING ERRORS RESISTANT PROGRAM FOR PREDICTING PROTEIN CODING REGIONS
#IN UNFINISHED cDNA. Journal of Bioinformatics and Computational Biology, 4(3):649 - 664, 2006.

sub angle_translation {
 my ($fasta_path, $output, $fasta_withorfs, $fasta_noorfs, $angle_dir, $angle_bit) = @_;
 my $temp_file = 'mfasta_fragment.temp';

 existe_arquivo($fasta_path);

 open(OUT, ">$output") or die "\nANGLE :: Error saving $output: $!\n";
 open(NOORF, ">$fasta_noorfs") or die "\nANGLE :: Error NoOrfs.fas $fasta_noorfs: $!\n";
 open(WITHORF, ">$fasta_withorfs") or die "\nANGLE :: Error WithOrfs.fas $fasta_withorfs: $!\n";
 select OUT; $|=1; select NOORF; $|=1; select WITHORF; $|=1; select STDOUT;
 my $notfirsttime = '0';
 my $sequencia;
 my $header;
 my $predicted_protein1;
 my $predicted_protein2;
 my $DP_score1;
 my $DP_score2;
 my $seq_size;

open(FASTA, $fasta_path) or die "\n$!";

while (my $linha = <FASTA>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  elsif ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco
##################---HEAD---######################
  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
      chomp $header;
       #SENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       $sequencia = uc $sequencia;
       print TEMP "$header\n";
       print TEMP "$sequencia";
       $seq_size = length ($sequencia);
       if ($seq_size => 30) { #we established min length = 80, so condition will always be met
        ($predicted_protein1, $DP_score1) = runparse_angle($angle_dir, $angle_bit, $temp_file);
        close TEMP;
        unlink($temp_file);
        #ANTISENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       print TEMP "$header\n";
       $sequencia = reverse($sequencia);
       print TEMP "$sequencia";
       ($predicted_protein2, $DP_score2) = runparse_angle($angle_dir, $angle_bit, $temp_file);
       close TEMP;
       unlink($temp_file);
       $sequencia = reverse($sequencia); #we need to undo the reversal to the seq
       
       if (($DP_score1 eq 'none') && ($DP_score2 eq 'none')) {print NOORF "$header\n"; print NOORF "$sequencia\n";}
       elsif (($DP_score1 ne 'none') && ($DP_score2 eq 'none')) {print OUT "$header\n"; print OUT "$predicted_protein1\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       elsif (($DP_score1 eq 'none') && ($DP_score2 ne 'none')) {print OUT "$header\n"; print OUT "$predicted_protein2\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       else {
        if ($DP_score1 => $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein1\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
        elsif ($DP_score1 < $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein2\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       }
       } #end if => 30
       else {print NOORF "$header\n"; print NOORF "$sequencia\n";}
  $sequencia = '';
  $header = $linha;
  }
##################---HEAD---######################
  elsif(($linha =~ /^\>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $header = $linha;
    $notfirsttime = '1';
  }
  elsif($linha !~ /^\>/) { #poderia ser else
   $sequencia .= $linha;
   chomp $sequencia;
   $sequencia =~ s/\s//g;
   $sequencia =~ s/\n//g;
  }

} #final while; processando ultima seq
##################---HEAD---######################
      chomp $header;
       #SENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       $sequencia = uc $sequencia;
       print TEMP "$header\n";
       print TEMP "$sequencia";
        $seq_size = length ($sequencia);
       if ($seq_size => 30) {
        ($predicted_protein1, $DP_score1) = runparse_angle($angle_dir, $angle_bit, $temp_file);
        close TEMP;
        unlink($temp_file);
        #ANTISENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       print TEMP "$header\n";
       $sequencia = reverse($sequencia);
       print TEMP "$sequencia";
       ($predicted_protein2, $DP_score2) = runparse_angle($angle_dir, $angle_bit, $temp_file);
       close TEMP;
       unlink($temp_file);
       $sequencia = reverse($sequencia);
       
       if (($DP_score1 eq 'none') && ($DP_score2 eq 'none')) {print NOORF "$header\n"; print NOORF "$sequencia\n";}
       elsif (($DP_score1 ne 'none') && ($DP_score2 eq 'none')) {print OUT "$header\n"; print OUT "$predicted_protein1\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       elsif (($DP_score1 eq 'none') && ($DP_score2 ne 'none')) {print OUT "$header\n"; print OUT "$predicted_protein2\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       else {
        if ($DP_score1 => $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein1\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
        elsif ($DP_score1 < $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein2\n"; print WITHORF "$header\n"; print WITHORF "$sequencia\n";}
       }
       } #end if => 30
       else {print NOORF "$header\n"; print NOORF "$sequencia\n";}
close FASTA;
close OUT;
close NOORF;
close WITHORF;
return;
} #end sub angle_translation
##################################

sub runparse_angle {
 my ($angle_dir, $angle_bit, $infile) = @_;
 my @angle_output = `$angle_dir/ANGLE-linux${angle_bit}DP -d $angle_dir/param-human -i $infile 2>&1`;

my %standard_genetic_code = (
 'TTT' => 'F',
 'TTC' => 'F',
 'TTA' => 'L',
 'TTG' => 'L',

 'CTT' => 'L',
 'CTC' => 'L',
 'CTA' => 'L',
 'CTG' => 'L',

 'ATT' => 'I',
 'ATC' => 'I',
 'ATA' => 'I',
 'ATG' => 'M',

 'GTT' => 'V',
 'GTC' => 'V',
 'GTA' => 'V',
 'GTG' => 'V',

 'TCT' => 'S',
 'TCC' => 'S',
 'TCA' => 'S',
 'TCG' => 'S',

 'CCT' => 'P',
 'CCC' => 'P',
 'CCA' => 'P',
 'CCG' => 'P',

 'ACT' => 'T',
 'ACC' => 'T',
 'ACA' => 'T',
 'ACG' => 'T',

 'GCT' => 'A',
 'GCC' => 'A',
 'GCA' => 'A',
 'GCG' => 'A',

 'TAT' => 'Y',
 'TAC' => 'Y',
 'TAA' => 'STOP',
 'TAG' => 'STOP',

 'CAT' => 'H',
 'CAC' => 'H',
 'CAA' => 'Q',
 'CAG' => 'Q',

 'AAT' => 'N',
 'AAC' => 'N',
 'AAA' => 'K',
 'AAG' => 'K',

 'GAT' => 'D',
 'GAC' => 'D',
 'GAA' => 'E',
 'GAG' => 'E',

 'TGT' => 'C',
 'TGC' => 'C',
 'TGA' => 'STOP',
 'TGG' => 'W',

 'CGT' => 'R',
 'CGC' => 'R',
 'CGA' => 'R',
 'CGG' => 'R',

 'AGT' => 'S',
 'AGC' => 'S',
 'AGA' => 'R',
 'AGG' => 'R',

 'GGT' => 'G',
 'GGC' => 'G',
 'GGA' => 'G',
 'GGG' => 'G',

);

my $get_seq_flag = 0;
my $starting_seq_flag = 0;
my $sequencia = 'undef';
my $orf_seq;
my $start;
my $stop;
my $DP_score;
ANGLEOUT:foreach my $line (@angle_output) {
 if (($get_seq_flag) && ($starting_seq_flag == 0)) {
  if ($line =~ /^\s*$/) {$starting_seq_flag = 1; $sequencia = '';}
 }
 if (($get_seq_flag) && ($starting_seq_flag)) {
  chomp $line; $line =~ s/\s//; $line =~ s/\n//;
  $sequencia .= $line;
 }
 if ($line =~ /DP SCORE: (-?\d+.?\d*)/) {
  $DP_score = $1; next ANGLEOUT;
 }
 if ($line =~ /START: (\d*) STOP: (\d*)/) {
   $start = $1; $stop = $2;
   #NO CODING SEQ FOUND::
   if (($start == 0) && ($stop == 0)) {$orf_seq = 'no_orf!'; $DP_score = 'none'; last ANGLEOUT;}
   #CODING SEQ FOUND::
   else {
    $get_seq_flag = 1;
    next ANGLEOUT;
   }
 }
}
if ($sequencia ne 'undef') {
 $sequencia =~ s/\{//; $sequencia =~ s/\}//; $sequencia =~ s/\(//; $sequencia =~ s/\)//; $sequencia =~ s/X//;
 $start = $start - 1;
  TRANSLATE:for($start; $start < $stop-2; $start += 3) {
   my $codon = substr($sequencia,$start,3);
   $codon = $standard_genetic_code{$codon};
   if ($codon eq 'STOP') {$codon = '';}
   $orf_seq .= $codon;
  }
}
return ($orf_seq, $DP_score);
} #end sub runparse_angle
######################################
sub checkFastaNT {
 my ($input_file) = @_;
 my $qtde_seqs = 0;
 my $flag_header = 0;
 my $flag_sequence = 0;
 my $invalid_chars = 0;

open(INPUT, "$input_file") or die "\nError opening $input_file: $!\n";

FILELOOP:while (my $line = <INPUT>) {

 if ($line =~ /^>/) {
  $qtde_seqs++;
  unless ($line =~ /^>\S+.+\n$/) {die "\nERROR: $input_file: Invalid FASTA file, violating header no. $qtde_seqs\n";} #sane header format
  if ($flag_header == 1) {die "\nERROR: $input_file : Invalid FASTA file, violating header no. $qtde_seqs\n";}
  $flag_header = 1;
  next FILELOOP;
 }
 elsif ($flag_header == 1) {
  if ($line =~ /^>/) {die "\nERROR: $input_file : Invalid FASTA file, violating header no. $qtde_seqs\n";}
  else {
   chomp $line;
   $flag_header = 0;
   $invalid_chars = ($line =~ tr/autgcnAUTGCN//c); if ($invalid_chars > 0) {die "\nERROR: $invalid_chars Invalid characters in $input_file header no. $qtde_seqs\n";} #remove unrecognized chars
   next FILELOOP;
  }
 }

} #end reading input
close INPUT;
return "$qtde_seqs";
} #sub CheckFastaNT
######################################
sub checkFastaPT {
 my ($input_file) = @_;
my $qtde_seqs = 0;
my $flag_header = 0;
my $flag_sequence = 0;
my $invalid_chars = 0;

open(INPUT, "$input_file") or die "\nError opening $input_file: $!\n";

FILELOOP:while (my $line = <INPUT>) {

 if ($line =~ /^>/) {
  $qtde_seqs++;
  unless ($line =~ /^>\S+.+\n$/) {die "\nERROR: $input_file : Invalid FASTA file, violating header no. $qtde_seqs\n";} #sane header format
  if ($flag_header == 1) {die "\nERROR: $input_file : Invalid FASTA file, violating header no. $qtde_seqs\n";}
  $flag_header = 1;
  next FILELOOP;
 }
 elsif ($flag_header == 1) {
  if ($line =~ /^>/) {die "\nERROR: $input_file : Invalid FASTA file, violating header no. $qtde_seqs\n";}
  else {
   chomp $line;
   $flag_header = 0;
   $invalid_chars = ($line =~ tr/IVLFCMAGTWSYPHEQDNKRXivlfcmagtwsypheqdnkrx//c); if ($invalid_chars > 0) {die "\nERROR: $invalid_chars Invalid characters in $input_file in header no. $qtde_seqs\n";} #remove unrecognized chars
   next FILELOOP;
  }
 }

} #end reading input
return "$qtde_seqs";

} #end sub checkFastaPT
######################################
sub fix_input {
 my ($input_file, $log_file) = @_;
 my $temp = "${input_file}_temp.fix";
 open(IN_FIX, $input_file) or die "\nFix input: $!\n";
 open(OUT_TEMP, ">$temp") or die "\nFix input: $!\n";
 select OUT_TEMP; $|=1; select LOGFILE; $|=1; select STDOUT;
 open(LOGFILE, ">>$log_file") or die "\nLog file: $!\n";
 my $invalid_chars = 0;
 my $notfirsttime = 0;
 my $sequencia;
 my $header;
 my $linha;
 my $error_count = 0;
 my $atleastoneseq = 0;
 my $count_N;

 while (my $linha = <IN_FIX>) {
  chomp $linha;
  #if ($linha =~ /#/) {next;} #ignorando partes comentadas#FRANKCODE
  if ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco#FRANKCODE
##################---HEAD---######################
  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
      chomp $header;
      $sequencia = uc $sequencia;
      $invalid_chars = ($sequencia =~ tr/AUTGCN//c);
      $count_N = ($sequencia =~ tr/N//);
      if ($invalid_chars > 0) {$error_count++; print LOGFILE "$header :: sequence ignored - contains $invalid_chars Invalid characters\n";}
      elsif (length($sequencia) < 80) {$error_count++; print LOGFILE "$header :: sequence ignored - smaller than 80 letters\n";}
      elsif (0.2 < ($count_N / length($sequencia))) {$error_count++; print LOGFILE "$header :: sequence ignored - too many \'N\' chars\n";}
      elsif ($header !~ /^>\S+.+/) {$error_count++; print LOGFILE "$header :: sequence ignored, invalid header.\n";} #sane header format
      else {                              #this one passed all FASTA quality tests
       $sequencia =~ tr/U/T/;             #reverse transcriptase
       $atleastoneseq = 1;
       print OUT_TEMP "$header\n";
       print OUT_TEMP "$sequencia\n";
      }
  $sequencia = '';
  $header = $linha;
  }
##################---HEAD---######################
  elsif(($linha =~ /^\>/) && ($notfirsttime == '0')) { #Loop pela primeira vez, apenas ler header
    $header = $linha;
    $notfirsttime = '1';
  }
  elsif($linha !~ /^\>/) { #poderia ser else
   $sequencia .= $linha;
   chomp $sequencia;
   $sequencia =~ s/\s//g;
   $sequencia =~ s/\n//g;
  }

} #final while; processando ultima seq
##################---HEAD---######################
      chomp $header;
      $sequencia = uc $sequencia;
      $invalid_chars = ($sequencia =~ tr/AUTGCN//c);
      $count_N = ($sequencia =~ tr/N//);
      if ($invalid_chars > 0) {$error_count++; print LOGFILE "$header :: sequence ignored - contains $invalid_chars Invalid characters\n";}
      elsif (length($sequencia) < 80) {$error_count++; print LOGFILE "$header :: sequence ignored - smaller than 80 letters\n";}
      elsif (0.2 < ($count_N / length($sequencia))) {$error_count++; print LOGFILE "$header :: sequence ignored - too many \'N\' chars\n";}
      elsif ($header !~ /^>\S+.+/) {$error_count++; print LOGFILE "$header :: sequence ignored, invalid header.\n";} #sane header format
      else {                              #this one passed all FASTA quality tests
       $sequencia =~ tr/U/T/;             #reverse transcriptase
       $atleastoneseq = 1;
       print OUT_TEMP "$header\n";
       print OUT_TEMP "$sequencia\n";
      }
if ($error_count > 0) {print LOGFILE "\nErrors: $error_count\n";} else {unlink $log_file;}
close IN_FIX;
close OUT_TEMP;
close LOGFILE;
 unlink ($input_file);
 rename($temp, $input_file);
unless ($atleastoneseq) {die "\nEXITING: all input sequences considered invalid. Please read $log_file for details.\n";}
return();
} #end sub fix_input
######################################
sub testAngle {
 my ($angle_dir) = @_;
 my @angle_output = `$angle_dir/ANGLE-linux32DP -d $angle_dir/param-human -i $angle_dir/sample.txt 2>&1`;
 if ($angle_output[2] =~ /FASTA HEADER/ ) {return '32';}
 else { #Let´s try the 64-bit flavour, then
  @angle_output = ();
  @angle_output = `$angle_dir/ANGLE-linux64DP -d $angle_dir/param-human -i $angle_dir/sample.txt 2>&1`;
  if ($angle_output[2] =~ /FASTA HEADER/ ) {return '64';}
  else {die "\nFATAL :: Unknown error in ANGLE program (32 and 64 bit versions). Please try installing PORTRAIT in another machine. ANGLE dump:\n\n @angle_output\n";}
 }
} #end sub testAngle
######################################
sub testAngleSecond {
 my ($angle_dir, $angle_bit) = @_;
 my @angle_output = `$angle_dir/ANGLE-linux${angle_bit}DP -d $angle_dir/param-human -i $angle_dir/sample.txt 2>&1`;
 if ($angle_output[2] =~ /FASTA HEADER/ ) {return;}
 else {die "\nFATAL :: Unknown error in ANGLE program ($angle_bit bit version). Please try installing PORTRAIT in another machine. ANGLE dump:\n\n @angle_output\n";}
}
######################################
sub testCast {
 my ($cast_dir) = @_;
 my @cast_output = `$cast_dir/cast -h 2>&1`;
 if ($cast_output[1] =~ /Compositional Bias Filtering Algorithm/) {return;}
 else {die "\nFATAL :: Unkwnown error in CAST program. If possible, please try installing PORTRAIT in another machine. CAST dump: \n\n @cast_output\n";}
} #end sub testCast
######################################
sub testLibsvm {
 my ($libsvm_dir) = @_;
 my @libsvm_output = `$libsvm_dir/svm-predict`;
 if ($libsvm_output[0] =~ /Usage: svm-predict/) {return;}
 else {die "\nFATAL :: Unkwnown error in LIBSVM program. If possible, please try reinstalling LIBSVM or installing PORTRAIT in another machine. LIBSVM dump: \n\n @libsvm_output\n";}
} #end sub testLibsvm
######################################


