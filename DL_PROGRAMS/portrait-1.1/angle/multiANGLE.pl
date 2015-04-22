#!/usr/bin/perl

use strict;

# RTARRIAL mar/2007

############################
# MULTIFASTA READER / FILTER
############################
# STATUS: FUNCIONAL


my $usage = "<path multifasta> <proteins output> <without orfs output>";
if ($#ARGV < 0) { die "\nUso: $0 $usage\n"; }

my $fasta_path = $ARGV[0];
my $output = $ARGV[1];
my $fasta_noorfs = $ARGV[2];
my $temp_file = 'mfasta_fragment.temp';
my $angle_dir = '/home/rtarrial/angle';

existe_arquivo($fasta_path);
nao_existe_arquivo($output);
nao_existe_arquivo($temp_file);
nao_existe_arquivo($fasta_noorfs);

open(OUT, ">$output") or die "\nErro ao salvar para $output: $!\n";
open(NOORF, ">$fasta_noorfs") or die "\nError NoOrfs.fas $fasta_noorfs: $!\n";
select OUT; $|=1; select NOORF; $|=1; select STDOUT;
my $notfirsttime = '0';
my $sequencia;
my $header;
my $predicted_protein1;
my $predicted_protein2;
my $DP_score1;
my $DP_score2;

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
       print TEMP "$header\n";
       print TEMP "$sequencia";
       ($predicted_protein1, $DP_score1) = runparse_angle($angle_dir, $temp_file);
       close TEMP;
       unlink($temp_file);
       #ANTISENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       print TEMP "$header\n";
       $sequencia = reverse($sequencia);
       print TEMP "$sequencia";
       ($predicted_protein2, $DP_score2) = runparse_angle($angle_dir, $temp_file);
       close TEMP;
       unlink($temp_file);

  if (($DP_score1 eq 'none') && ($DP_score2 eq 'none')) {print NOORF "$header\n"; print NOORF "$sequencia\n";}
  elsif (($DP_score1 ne 'none') && ($DP_score2 eq 'none')) {print OUT "$header\n"; print OUT "$predicted_protein1\n";}
  elsif (($DP_score1 eq 'none') && ($DP_score2 ne 'none')) {print OUT "$header\n"; print OUT "$predicted_protein2\n";}
  else {
   if ($DP_score1 => $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein1\n";}
   elsif ($DP_score1 < $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein2\n";}
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
       #SENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       print TEMP "$header\n";
       print TEMP "$sequencia";
       ($predicted_protein1, $DP_score1) = runparse_angle($angle_dir, $temp_file);
       close TEMP;
       unlink($temp_file);
       #ANTISENSE PREDICTION
       open(TEMP, ">$temp_file") or die "\nError saving to temp file $temp_file: $!\n";
       select TEMP; $|=1; select STDOUT;
       print TEMP "$header\n";
       $sequencia = reverse($sequencia);
       print TEMP "$sequencia";
       ($predicted_protein2, $DP_score2) = runparse_angle($angle_dir, $temp_file);
       close TEMP;
       unlink($temp_file);
  #BACK TO SENSE
  $sequencia = reverse($sequencia);
  #SEQ WITH BIGGEST DP WINS!
  if (($DP_score1 eq 'none') && ($DP_score2 eq 'none')) {print NOORF "$header\n"; print NOORF "$sequencia\n";}
  elsif (($DP_score1 ne 'none') && ($DP_score2 eq 'none')) {print OUT "$header\n"; print OUT "$predicted_protein1\n";}
  elsif (($DP_score1 eq 'none') && ($DP_score2 ne 'none')) {print OUT "$header\n"; print OUT "$predicted_protein2\n";}
  else {
   if ($DP_score1 => $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein1\n";}
   elsif ($DP_score1 < $DP_score2) {print OUT "$header\n"; print OUT "$predicted_protein2\n";}
  }
close FASTA;
close OUT;
close NOORF;
print "\nDone. Proteins saved in $output, FASTA without ORFs saved in $fasta_noorfs\n";
exit;


##################################
sub existe_arquivo {
        my ($arquivo) = @_;
        if (-e $arquivo) { ; } else { die "\nArquivo $arquivo nao existe!\n\n"; }
}

##################################

sub nao_existe_arquivo {
        my ($arquivo) = @_;
        if (-e $arquivo) { die "\nArquivo de saida $arquivo ja existe!\n";} else { ; }
}

##################################

sub abre_arquivo {
        my ($arquivo) = @_;
        my @conteudo = ();
        open (ARQ, $arquivo) or die "\nNao pude abrir $arquivo.\n";
        @conteudo = <ARQ>;
        close ARQ;
        return @conteudo;
}

##################################
#------------------------------------------------------
#Now Calculating...

#FASTA HEADER: seq1
#SEQ LENGTH: 1027
#GC PRESSURE: 48.9%
#DP SCORE: 289.3
#START: 183 STOP: 1023
#INSERTION_ERR: A@1023

#AA{ATGTAGACAACAAGCGGAAAGAAGAGGCAGTCCAATTCCTCCAACAATGCTGCCGCC
#CTGCTGAACTCCGCCTGTCTGGAACTCTTCGAATTCGTCAAACGCGAAAATATTAAGGCG
#ATA(A)}TCAT
#------------------------------------------------------
#Now Calculating...

#FASTA HEADER: seq2
#SEQ LENGTH: 66
#GC PRESSURE: 60.6%
#START: 0 STOP: 0

#NO CODING REGION WAS FOUND IN THIS SEQUENCE.



sub runparse_angle {
 my ($angle_dir, $infile) = @_;
 if ($angle_dir =~ /(.*)\/$/) {$angle_dir = $1;} #getting rid of last slash
 my @angle_output = `$angle_dir/ANGLE-linux32DP -i $infile 2>&1`;

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
