#!usr/bin/perl
#------------------------------------------------------
#Now Calculating...

#FASTA HEADER: seq1
#SEQ LENGTH: 1027
#GC PRESSURE: 48.9%
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

my $angle_dir = '/home/rtarrial/angle';
my @angle_output = `$angle_dir/ANGLE-linux32 -i fragment.fas 2>&1`;

my $get_seq_flag = 0;
my $starting_seq_flag = 0;
my $sequencia = 'undef';
my $orf_seq;
my $start;
my $stop;
ANGLEOUT:foreach my $line (@angle_output) {
 if (($get_seq_flag) && ($starting_seq_flag == 0)) {
  if ($line =~ /^\s*$/) {$starting_seq_flag = 1; $sequencia = '';}
 }
 if (($get_seq_flag) && ($starting_seq_flag)) {
  chomp $line; $line =~ s/\s//; $line =~ s/\n//;
  $sequencia .= $line;
 }
 if ($line =~ /START: (\d*) STOP: (\d*)/) {
   $start = $1; $stop = $2;
   #NO CODING SEQ FOUND::
   if (($start == 0) && ($stop == 0)) {$orf_seq = 'no_orf!'; last ANGLEOUT;}
   #CODING SEQ FOUND::
   else {
    $get_seq_flag = 1;
   }
 }
}
$sequencia =~ s/\{//;$sequencia =~ s/\}//;$sequencia =~ s/\(//;$sequencia =~ s/\)//;$sequencia =~ s/X//;
my $orf_seq;
$start = $start - 1;
if ($sequencia ne 'undef') {
  TRANSLATE:for($start; $start < $stop-2; $start += 3) {
   my $codon = substr($sequencia,$start,3);
   $codon = $standard_genetic_code{$codon};
   if ($codon eq 'STOP') {$codon = '';}
   $orf_seq .= $codon;
  }
}
print "PTN: $orf_seq\n";

