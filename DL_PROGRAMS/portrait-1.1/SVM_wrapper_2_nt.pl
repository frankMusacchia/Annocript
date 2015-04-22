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
#   PORTRAIT Prediction of transcriptomic ncRNA by ab initio methods      #                                                        
#                                                                         #
#                                                                         #
#   Description:                                                          #
#      an algorithm for identifying non-coding RNAs on a transcriptome    #
#      context.                                                           #
#                                                                         #
#   Developed by:                                                         #
#                                                                         #
###########################################################################

use warnings;
use strict;

$| = 1; #AUTOFLUSH

#RTARRIAL 2006/2007
#SPARSE Format

if ($#ARGV < 0) {
 print "\nThis script is called by SVM_inputter_auto.pl. Please run it instead.\n";
 exit;
}

my $indir = $ARGV[0];
my $classe = $ARGV[1];

if (defined($classe)) {
 if ($classe == '1') { ; }
 elsif ($classe == '-1') { ; }
 elsif ($classe == '0') { ; }
 else {die "\nClass options: \n1 = coding, \n-1 = non-coding, \n0 = test set.\n";}
} #end defined
else {die "\nClass options: \n1 = coding, \n-1 = non-coding, \n0 = test set.\n";}

my $fastafile = "${indir}";

my $cria_rnasize = 0;
my $cria_nt1 = 0;
my $cria_nt2 = 0;
my $cria_nt3 = 0;
my $rnasizefile;
my $nt1file;
my $nt2file;
my $nt3file;

existe_arquivo($fastafile);
nao_existe_arquivo("${indir}.SVM");

if (-e "${indir}_N.rnasize") { $rnasizefile = "${indir}_N.rnasize"; } else { $cria_rnasize  = 1; }
if (-e "${indir}_N.nt1") { $nt1file = "${indir}_N.nt1"; } else { $cria_nt1 = 1; }
if (-e "${indir}_N.nt2") { $nt2file = "${indir}_N.nt2"; } else { $cria_nt2 = 1; }
if (-e "${indir}_N.nt3") { $nt3file = "${indir}_N.nt3"; } else { $cria_nt3 = 1; }

# Todos os arquivos existem, agora abrir um a um e parsear

#####################
# PROCESSING
#####################

#COUNT AMOUNT OF SEQS IN FASTAFILE
my @fastafile = abre_arquivo($fastafile);
my $qtde_seqs = 0;
foreach my $linha (@fastafile) {
 if ($linha =~ /^>/) { $qtde_seqs++; }
}
#-=-=-=-=-#
mfasta_reader_nucleotide($indir, $fastafile, $cria_rnasize, $cria_nt1, $cria_nt2, $cria_nt3);
$rnasizefile = "${indir}_N.rnasize";
$nt1file = "${indir}_N.nt1";
$nt2file = "${indir}_N.nt2";
$nt3file = "${indir}_N.nt3";
#
my @tamanho_transcrito = abre_arquivo($rnasizefile);
my $amount_tamanho = @tamanho_transcrito;
if ($amount_tamanho != $qtde_seqs) { die "\n Fatal Error: Amount of transcript sizes ($amount_tamanho) differ from FASTA ($qtde_seqs)\n"; }
#-=-=-=-=-#
my @nucleotideos = abre_arquivo($nt1file);
my $amount_nt = @nucleotideos;
if ($amount_nt != $qtde_seqs) { die "\n Fatal Error: Amount of nucleotides ($amount_nt) differ from FASTA ($qtde_seqs)\n"; }
#-=-=-=-=-#
my @dinucleotideos = abre_arquivo($nt2file);
my $amount_nt2 = @dinucleotideos;
if ($amount_nt2 != $qtde_seqs) { die "\n Fatal Error: Amount of dinucleotides ($amount_nt2) differ from FASTA ($qtde_seqs)\n"; }
#-=-=-=-=-#
my @trinucleotideos = abre_arquivo($nt3file);
my $amount_nt3 = @trinucleotideos;
if ($amount_nt3 != $qtde_seqs) { die "\n Fatal Error: Amount of trinucleotides ($amount_nt3) differ from FASTA ($qtde_seqs)\n"; }

###########
###############################################
# PROCESSANDO OUTPUT
################################
###############################################

open(OUT, ">${indir}.SVM") or die "\nError saving ${indir}.SVM: $! \n";

#ORDEM DOS ATRIBUTOS
#Nucleotideos: 1-4
#Dinucleotideos: 5-20
#Trinucleotideos: 21-84
#Aminoacidos: 85-104
#Tamanho ORF: 105-108
#CAST: 109
#SOAP: 110
#IEP: 111
#FoldIndex: 112
#PSSCAN: 113
#BLAST: 114
#SSpro/ACCpro: 115-118


for (my $j='0'; $j <= ($qtde_seqs-1); $j++) {
 my $firstone;
 my $seqatrib;

 #nt1
 my $nt1;
 $nt1 = $nucleotideos[$j];
 chomp $nt1;
 my @nt1_svm = ();
 @nt1_svm = split(/ /, $nt1);
 $nt1 = '';
#SPARSE formatting
 $firstone = 1;
 for (my $seqatrib = 0; $seqatrib <= 3; $seqatrib++){
  if ($firstone) {
   if ($nt1_svm[($seqatrib)] > 0){
    $nt1 .= $seqatrib+1 . ":$nt1_svm[$seqatrib]";
    $firstone = 0;
   }
  }
  elsif ($nt1_svm[($seqatrib)] > 0){
    $nt1 .= " ";  #nao consegui incorporar " " na string abaixo!
    $nt1 .= $seqatrib+1 . ":$nt1_svm[$seqatrib]";
  }
 } #end for
#$nt1 = "1:$nt1_svm[0] 2:$nt1_svm[1] 3:$nt1_svm[2] 4:$nt1_svm[3]";
 
  #nt2
 my $nt2;
 $nt2 = $dinucleotideos[$j];
 chomp $nt2;
 my @nt2_svm = ();
 @nt2_svm = split(/ /, $nt2);

#SPARSE formatting
 $nt2 = '';
 $firstone = 1;
 for (my $seqatrib = 0; $seqatrib <= 15; $seqatrib++){
  if ($firstone) {
   if ($nt2_svm[($seqatrib)] > 0){
    $nt2 .= $seqatrib+5 . ":$nt2_svm[$seqatrib]";
    $firstone = 0;
   }
  }
  elsif ($nt2_svm[($seqatrib)] > 0){
    $nt2 .= " ";  #nao consegui incorporar " " na string abaixo!
    $nt2 .= $seqatrib+5 . ":$nt2_svm[$seqatrib]";
  }
 }#end for
#$nt2 = "5:$nt2_svm[0] 6:$nt2_svm[1] 7:$nt2_svm[2] 8:$nt2_svm[3] 9:$nt2_svm[4] 10:$nt2_svm[5] 11:$nt2_svm[6] 12:$nt2_svm[7] 13:$nt2_svm[8] 14:$nt2_svm[9] 15:$nt2_svm[10] 16:$nt2_svm[11] 17:$nt2_svm[12] 18:$nt2_svm[13] 19:$nt2_svm[14] 20:$nt2_svm[15]";

  #nt3
 my $nt3;
 $nt3 = $trinucleotideos[$j];
 chomp $nt3;
 my @nt3_svm = ();
 @nt3_svm = split(/ /, $nt3);
#SPARSE formatting
 $nt3 = '';
 $firstone = 1;
 for (my $seqatrib = 0; $seqatrib <= 63; $seqatrib++){
  if ($firstone) {
   if ($nt3_svm[($seqatrib)] > 0){
    $nt3 .= $seqatrib+21 . ":$nt3_svm[$seqatrib]";
    $firstone = 0;
   }
  }
  elsif ($nt3_svm[($seqatrib)] > 0){
    $nt3 .= " ";  #nao consegui incorporar " " na string abaixo!
    $nt3 .= $seqatrib+21 . ":$nt3_svm[$seqatrib]";
  }
 }#end for

#$nt3 = "21:$nt3_svm[0] 22:$nt3_svm[1] 23:$nt3_svm[2] 24:$nt3_svm[3] 25:$nt3_svm[4] 26:$nt3_svm[5] 27:$nt3_svm[6] 28:$nt3_svm[7] 29:$nt3_svm[8] 30:$nt3_svm[9] 31:$nt3_svm[10] 32:$nt3_svm[11] 33:$nt3_svm[12] 34:$nt3_svm[13] 35:$nt3_svm[14] 36:$nt3_svm[15] 37:$nt3_svm[16] 38:$nt3_svm[17] 39:$nt3_svm[18] 40:$nt3_svm[19] 41:$nt3_svm[20] 42:$nt3_svm[21] 43:$nt3_svm[22] 44:$nt3_svm[23] 45:$nt3_svm[24] 46:$nt3_svm[25] 47:$nt3_svm[26] 48:$nt3_svm[27] 49:$nt3_svm[28] 50:$nt3_svm[29] 51:$nt3_svm[30] 52:$nt3_svm[31] 53:$nt3_svm[32] 54:$nt3_svm[33] 55:$nt3_svm[34] 56:$nt3_svm[35] 57:$nt3_svm[36] 58:$nt3_svm[37] 59:$nt3_svm[38] 60:$nt3_svm[39] 61:$nt3_svm[40] 62:$nt3_svm[41] 63:$nt3_svm[42] 64:$nt3_svm[43] 65:$nt3_svm[44] 66:$nt3_svm[45] 67:$nt3_svm[46] 68:$nt3_svm[47] 69:$nt3_svm[48] 70:$nt3_svm[49] 71:$nt3_svm[50] 72:$nt3_svm[51] 73:$nt3_svm[52] 74:$nt3_svm[53] 75:$nt3_svm[54] 76:$nt3_svm[55] 77:$nt3_svm[56] 78:$nt3_svm[57] 79:$nt3_svm[58] 80:$nt3_svm[59] 81:$nt3_svm[60] 82:$nt3_svm[61] 83:$nt3_svm[62] 84:$nt3_svm[63]";

 #TAMANHO RNA
 my $tamanho_rna_svm;
 $tamanho_rna_svm = categoriza_rnasize($tamanho_transcrito[$j]);



print OUT "$classe $nt1 $nt2 $nt3 $tamanho_rna_svm\n";
}
close OUT;
checa_libsvm("${indir}.SVM");

##################################
###############################
#      SUBROTINAS             #
###############################
##################################
sub existe_arquivo {
        my ($arquivo) = @_;
        if (-e $arquivo) { ; } else { die "\nCould not find file $arquivo.\n\n"; }
}

##################################

sub nao_existe_arquivo {
        my ($arquivo) = @_;
        if (-e $arquivo) { die "\nFile $arquivo already exists.\n";} else { ; }
}

##################################

sub abre_arquivo {
        my ($arquivo) = @_;
        my @conteudo = ();
        open (ARQ, $arquivo) or die "\n(sub abre_arquivo)::Could not open $arquivo.\n";
        select ARQ; $| = 1; select STDOUT;
        @conteudo = <ARQ>;
        close ARQ;
        return @conteudo;
}
###############################
#########################
#Checa formato libSVM
#########################
# Executa diversos testes para ver se o arquivo gerado pelo SVM_inputter
# está no formato exigido pelo libSVM.
# Espera formato SPARSE.
#1:0.25222882615156 2:0.256686478454681 3:0.24739970282318 4:0.243684992570579 5:0.0705794947994056 6:0.0676077265973254
#1:0.26115760111576 2:0.234309623430962 3:0.239888423988842 4:0.264644351464435 5:0.0801952580195258 6:0.0550906555090655
#[-11nothing] integer:1>real>0 integer:1>real>0 ...
sub checa_libsvm {
my ($input) = @_;

open(IN, $input) or die "\n$input -- $!\n";
my $checa_tipo = 1;
my $tipo = 'nulo';
my $instances = 0;
my $atributos = 0;

while (my $linha = <IN>) {
 chomp $linha;
 $instances++;
 if ($checa_tipo) {
  $linha =~ /^-?\d?:?/;
  if ($& =~ /0/) {$tipo = 'teste';}
  elsif ($& == '-1') {$tipo = 'negativo de treinamento';}
  elsif ($& == '1') {$tipo = 'positivo de treinamento';}
  else {die "\nErro no arquivo: Classe do arquivo nao reconhecida: $&\n";}
  $checa_tipo = 0;
 }
#Manter quantidade constante de atributos entre instancias -- INVALIDO PARA SPARSE
# if (($linha =~ tr/://) != $atributos) {die "\nQuantidade de atributos diferentes na linha $instances.\n";}

while ($linha =~ /\D/g) {
 if (($& ne ":") && ($& ne ".") && ($& ne " ") && ($& ne "-")) {die "\nCaractere nao-numerico: \"$&\" na linha $instances\n";}
}

# my $atrib = 0;
 while ($linha =~ /(\d+):(\d+\.?\d*)/g) {
#NAO FUNCIONA EM FORMATO SPARSE :: if ($1 != $atrib) {die "\nNumero sequencial errado de atributo: esperava $atrib, obtive $1 na linha $instances.\n";}
  if ($1 =~ /\D+/g) {die "\nCaractere nao-numerico na indicacao de atributo: \"$&\" na linha $instances\n";}
  if (($2 > 0) && ($2 <= 1)) { ; } else {die "\nValor de atributo fora do range \[0,1\]: valor $2 na linha $instances.\n";}
  if ($2 =~ /\.\D+/g) {die "\nCaractere nao-numerico no valor de atributo: \"$&\" na linha $instances\n";}
 }
# if ($atrib != $atributos) {die "\nQuantidade de atributos diferente na linha $instances.\n";}

} #fim do arquivo
return;
} #end sub checa_libsvm
########################
# ONE MFASTA_READER TO EXTRACT ATTRIBUTES FROM THEM ALL!!
sub mfasta_reader_nucleotide {
 my ($indir, $fasta_path, $cria_rnasize, $cria_nt1, $cria_nt2, $cria_nt3) = @_;
 if ($cria_rnasize) { my $out_rnasize = "${fasta_path}_N.rnasize"; nao_existe_arquivo($out_rnasize); open(OUT_RNASIZE, ">$out_rnasize") or die "\nOUT_RNASIZE: $out_rnasize - $!\n"; select OUT_RNASIZE; $|=1;}
 if ($cria_nt1) { my $out_nt1 = "${fasta_path}_N.nt1"; nao_existe_arquivo($out_nt1); open(OUT_NT1, ">$out_nt1") or die "\nOUT_NT1: $out_nt1 - $!\n"; select OUT_NT1; $|=1;}
 if ($cria_nt2) { my $out_nt2 = "${fasta_path}_N.nt2"; nao_existe_arquivo($out_nt2); open(OUT_NT2, ">$out_nt2") or die "\nOUT_NT2: $out_nt2 - $!\n";select OUT_NT2; $|=1;}
 if ($cria_nt3) { my $out_nt3 = "${fasta_path}_N.nt3"; nao_existe_arquivo($out_nt3); open(OUT_NT3, ">$out_nt3") or die "\nOUT_NT3: $out_nt3 - $!\n";select OUT_NT3; $|=1;}
 select STDOUT;

 my $notfirsttime = '0';
 my $sequencia;
 my $header;
 my $tamanho_seq;
 my $temp_value;

 open(FASTA, $fasta_path) or die "\n$!\n";

 while (my $linha = <FASTA>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  elsif ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco
##################---HEAD---######################
  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
     if ($cria_rnasize) {$tamanho_seq = length($sequencia); print OUT_RNASIZE "$tamanho_seq\n";}
     if ($cria_nt1) {$temp_value = return_nt1($sequencia); print OUT_NT1 "$temp_value\n";}
     if ($cria_nt2) {$temp_value = return_nt2($sequencia); print OUT_NT2 "$temp_value\n";}
     if ($cria_nt3) {$temp_value = return_nt3($sequencia); print OUT_NT3 "$temp_value\n";}
  $temp_value = '';
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
     if ($cria_rnasize) {$tamanho_seq = length($sequencia); print OUT_RNASIZE "$tamanho_seq\n";}
     if ($cria_nt1) {$temp_value = return_nt1($sequencia); print OUT_NT1 "$temp_value\n";}
     if ($cria_nt2) {$temp_value = return_nt2($sequencia); print OUT_NT2 "$temp_value\n";}
     if ($cria_nt3) {$temp_value = return_nt3($sequencia); print OUT_NT3 "$temp_value\n";}
close FASTA;
return;
} #end sub mfasta_reader
#################################
sub return_nt1 {
my ($sequencia) = @_;
my $qtde_A = ($sequencia =~ tr/Aa//);
my $qtde_T = ($sequencia =~ tr/Tt//);
my $qtde_G = ($sequencia =~ tr/Gg//);
my $qtde_C = ($sequencia =~ tr/Cc//);
#$qtde_N = ($sequencia =~ tr/Nn//);
#$qtde_X = ($sequencia =~ tr/Xx//);
my $tot = $qtde_A+$qtde_T+$qtde_G+$qtde_C;
# comentar a linha seguinte faz com que o output deixe de ser normalizado
$qtde_A = sprintf("%.17f", $qtde_A/$tot); $qtde_T = sprintf("%.17f", $qtde_T/$tot); $qtde_G = sprintf("%.17f", $qtde_G/$tot); $qtde_C = sprintf("%.17f", $qtde_C/$tot);
return "$qtde_A $qtde_T $qtde_G $qtde_C";
} #end sub return_nt1
#################################
sub return_nt2 {
my ($sequencia) = @_;
my $codon;
my $tot;
my $tamanho_seq;

my $AA = '0'; my $AT = '0'; my $AG = '0'; my $AC = '0';
my $TA = '0'; my $TT = '0'; my $TG = '0'; my $TC = '0';
my $GA = '0'; my $GT = '0'; my $GG = '0'; my $GC = '0';
my $CA = '0'; my $CT = '0'; my $CG = '0'; my $CC = '0';
$tamanho_seq = length($sequencia);

     for(my $i=0; $i < $tamanho_seq; $i += 2) {
       $codon = substr($sequencia,$i,2);
       $codon = uc $codon; #converte pra maiuscula
       if ($sequencia =~ /^\s*$/) { next;} #ignorando espaços em branco
       if ($codon !~ /[ATGC][ATGC]/) {next;} #descartar codons com X, N...

       if ($codon =~ /AA/) { $AA++;}
       elsif ($codon =~ /AT/) { $AT++;}
       elsif ($codon =~ /AG/) { $AG++;}
       elsif ($codon =~ /AC/) { $AC++;}

       elsif ($codon =~ /TA/) { $TA++;}
       elsif ($codon =~ /TT/) { $TT++;}
       elsif ($codon =~ /TG/) { $TG++;}
       elsif ($codon =~ /TC/) { $TC++;}

       elsif ($codon =~ /GA/) { $GA++;}
       elsif ($codon =~ /GT/) { $GT++;}
       elsif ($codon =~ /GG/) { $GG++;}
       elsif ($codon =~ /GC/) { $GC++;}

       elsif ($codon =~ /CA/) { $CA++;}
       elsif ($codon =~ /CT/) { $CT++;}
       elsif ($codon =~ /CG/) { $CG++;}
       elsif ($codon =~ /CC/) { $CC++;}

      } #end for

     $tot = '0';
     $tot = $AA+$AT+$AG+$AC+$TA+$TT+$TG+$TC+$GA+$GT+$GG+$GC+$CA+$CT+$CG+$CC;
     $AA = sprintf("%.17f", $AA/$tot); $AT = sprintf("%.17f", $AT/$tot); $AG = sprintf("%.17f", $AG/$tot); $AC = sprintf("%.17f", $AC/$tot);
     $TA = sprintf("%.17f", $TA/$tot); $TT = sprintf("%.17f", $TT/$tot); $TG = sprintf("%.17f", $TG/$tot); $TC = sprintf("%.17f", $TC/$tot);
     $GA = sprintf("%.17f", $GA/$tot); $GT = sprintf("%.17f", $GT/$tot); $GG = sprintf("%.17f", $GG/$tot); $GC = sprintf("%.17f", $GC/$tot);
     $CA = sprintf("%.17f", $CA/$tot); $CT = sprintf("%.17f", $CT/$tot); $CG = sprintf("%.17f", $CG/$tot); $CC = sprintf("%.17f", $CC/$tot);

return "$AA $AT $AG $AC $TA $TT $TG $TC $GA $GT $GG $GC $CA $CT $CG $CC";

} #end sub return_nt2
###############################
sub return_nt3 {
my ($sequencia) = @_;
my $codon;
my $tot = '0';
my $erros = '0'; #debug purposes only

   my $AAA = '0'; my $AAT = '0'; my $AAG = '0'; my $AAC = '0';
   my $ATA = '0'; my $ATT = '0'; my $ATG = '0'; my $ATC = '0';
   my $AGA = '0'; my $AGT = '0'; my $AGG = '0'; my $AGC = '0';
   my $ACA = '0'; my $ACT = '0'; my $ACG = '0'; my $ACC = '0';

   my $TAA = '0'; my $TAT = '0'; my $TAG = '0'; my $TAC = '0';
   my $TTA = '0'; my $TTT = '0'; my $TTG = '0'; my $TTC = '0';
   my $TGA = '0'; my $TGT = '0'; my $TGG = '0'; my $TGC = '0';
   my $TCA = '0'; my $TCT = '0'; my $TCG = '0'; my $TCC = '0';

   my $GAA = '0'; my $GAT = '0'; my $GAG = '0'; my $GAC = '0';
   my $GTA = '0'; my $GTT = '0'; my $GTG = '0'; my $GTC = '0';
   my $GGA = '0'; my $GGT = '0'; my $GGG = '0'; my $GGC = '0';
   my $GCA = '0'; my $GCT = '0'; my $GCG = '0'; my $GCC = '0';

   my $CAA = '0'; my $CAT = '0'; my $CAG = '0'; my $CAC = '0';
   my $CTA = '0'; my $CTT = '0'; my $CTG = '0'; my $CTC = '0';
   my $CGA = '0'; my $CGT = '0'; my $CGG = '0'; my $CGC = '0';
   my $CCA = '0'; my $CCT = '0'; my $CCG = '0'; my $CCC = '0';

          for(my $i=0; $i < (length($sequencia) - 2) ; $i += 3) {

          $codon = substr($sequencia,$i,3);
          $codon = uc $codon; #converte pra maiuscula
          if ($codon !~ /[ATGC][ATGC][ATGC]/) {next;} #descartar codons com X, N...

          if ($codon =~ /A../) {
            if ($codon =~ /AA./) {
              if ($codon eq 'AAA') { $AAA++;}
              elsif ($codon eq 'AAT') { $AAT++;}
              elsif ($codon eq 'AAG') { $AAG++;}
              elsif ($codon eq 'AAC') { $AAC++;}
            }
            elsif ($codon =~ /AT./) {
              if ($codon eq 'ATA') { $ATA++;}
              elsif ($codon eq 'ATT') { $ATT++;}
              elsif ($codon eq 'ATG') { $ATG++;}
              elsif ($codon eq 'ATC') { $ATC++;}
            }
            elsif ($codon =~ /AG./) {
              if ($codon eq 'AGA') { $AGA++;}
              elsif ($codon eq 'AGT') { $AGT++;}
           elsif ($codon eq 'AGG') { $AGG++;}
        elsif ($codon eq 'AGC') { $AGC++;}
        }
            elsif ($codon =~ /AC./) {
              if ($codon eq 'ACA') { $ACA++;}
              elsif ($codon eq 'ACT') { $ACT++;}
              elsif ($codon eq 'ACG') { $ACG++;}
              elsif ($codon eq 'ACC') { $ACC++;}
            }
         }
###############
###############

         elsif ($codon =~ /T../) {
            if ($codon =~ /TA./) {
              if ($codon eq 'TAA') { $TAA++;}
              elsif ($codon eq 'TAT') { $TAT++;}
              elsif ($codon eq 'TAG') { $TAG++;}
              elsif ($codon eq 'TAC') { $TAC++;}
            }
            elsif ($codon =~ /TT./) {
              if ($codon eq 'TTA') { $TTA++;}
              elsif ($codon eq 'TTT') { $TTT++;}
              elsif ($codon eq 'TTG') { $TTG++;}
              elsif ($codon eq 'TTC') { $TTC++;}
            }
            elsif ($codon =~ /TG./) {
              if ($codon eq 'TGA') { $TGA++;}
              elsif ($codon eq 'TGT') { $TGT++;}
              elsif ($codon eq 'TGG') { $TGG++;}
              elsif ($codon eq 'TGC') { $TGC++;}
            }
            elsif ($codon =~ /TC./) {
              if ($codon eq 'TCA') { $TCA++;}
              elsif ($codon eq 'TCT') { $TCT++;}
              elsif ($codon eq 'TCG') { $TCG++;}
              elsif ($codon eq 'TCC') { $TCC++;}
            }

         }

        ###############
        ###############

         elsif ($codon =~ /G../) {
            if ($codon =~ /GA./) {
              if ($codon eq 'GAA') { $GAA++;}
              elsif ($codon eq 'GAT') { $GAT++;}
              elsif ($codon eq 'GAG') { $GAG++;}
              elsif ($codon eq 'GAC') { $GAC++;}
            }
            elsif ($codon =~ /GT./) {
              if ($codon eq 'GTA') { $GTA++;}
              elsif ($codon eq 'GTT') { $GTT++;}
              elsif ($codon eq 'GTG') { $GTG++;}
              elsif ($codon eq 'GTC') { $GTC++;}
            }
            elsif ($codon =~ /GG./) {
              if ($codon eq 'GGA') { $GGA++;}
              elsif ($codon eq 'GGT') { $GGT++;}
              elsif ($codon eq 'GGG') { $GGG++;}
              elsif ($codon eq 'GGC') { $GGC++;}
            }
            elsif ($codon =~ /GC./) {
              if ($codon eq 'GCA') { $GCA++;}
              elsif ($codon eq 'GCT') { $GCT++;}
              elsif ($codon eq 'GCG') { $GCG++;}
              elsif ($codon eq 'GCC') { $GCC++;}
            }
         }

###############
###############

         elsif ($codon =~ /C../) {
            if ($codon =~ /CA./) {
              if ($codon eq 'CAA') { $CAA++;}
              elsif ($codon eq 'CAT') { $CAT++;}
              elsif ($codon eq 'CAG') { $CAG++;}
              elsif ($codon eq 'CAC') { $CAC++;}
            }
            elsif ($codon =~ /CT./) {
              if ($codon eq 'CTA') { $CTA++;}
              elsif ($codon eq 'CTT') { $CTT++;}
              elsif ($codon eq 'CTG') { $CTG++;}
              elsif ($codon eq 'CTC') { $CTC++;}
            }
            elsif ($codon =~ /CG./) {
              if ($codon eq 'CGA') { $CGA++;}
              elsif ($codon eq 'CGT') { $CGT++;}
              elsif ($codon eq 'CGG') { $CGG++;}
              elsif ($codon eq 'CGC') { $CGC++;}
            }
            elsif ($codon =~ /CC./) {
              if ($codon eq 'CCA') { $CCA++;}
              elsif ($codon eq 'CCT') { $CCT++;}
              elsif ($codon eq 'CCG') { $CCG++;}
              elsif ($codon eq 'CCC') { $CCC++;}
            }

         }
         else { $erros++; }

        } #end for loop (individual)

     $tot = '0';
     $tot = $AAA+$AAT+$AAG+$AAC+$ATA+$ATT+$ATG+$ATC+$AGA+$AGT+$AGG+$AGC+$ACA+$ACT+$ACG+$ACC+$TAA+$TAT+$TAG+$TAC+$TTA+$TTT+$TTG+$TTC+$TGA+$TGT+$TGG+$TGC+$TCA+$TCT+$TCG+$TCC+$GAA+$GAT+$GAG+$GAC+$GTA+$GTT+$GTG+$GTC+$GGA+$GGT+$GGG+$GGC+$GCA+$GCT+$GCG+$GCC+$CAA+$CAT+$CAG+$CAC+$CTA+$CTT+$CTG+$CTC+$CGA+$CGT+$CGG+$CGC+$CCA+$CCT+$CCG+$CCC;

  #Se as seguintes linhas forem comentadas, o output deixa de ser normalizado:

        $AAA = sprintf("%.17f", $AAA/$tot); $AAT = sprintf("%.17f", $AAT/$tot); $AAG = sprintf("%.17f", $AAG/$tot); $AAC = sprintf("%.17f", $AAC/$tot);
        $ATA = sprintf("%.17f", $ATA/$tot); $ATT = sprintf("%.17f", $ATT/$tot); $ATG = sprintf("%.17f", $ATG/$tot); $ATC = sprintf("%.17f", $ATC/$tot);
        $AGA = sprintf("%.17f", $AGA/$tot); $AGT = sprintf("%.17f", $AGT/$tot); $AGG = sprintf("%.17f", $AGG/$tot); $AGC = sprintf("%.17f", $AGC/$tot);
        $ACA = sprintf("%.17f", $ACA/$tot); $ACT = sprintf("%.17f", $ACT/$tot); $ACG = sprintf("%.17f", $ACG/$tot); $ACC = sprintf("%.17f", $ACC/$tot);

        $TAA = sprintf("%.17f", $TAA/$tot); $TAT = sprintf("%.17f", $TAT/$tot); $TAG = sprintf("%.17f", $TAG/$tot); $TAC = sprintf("%.17f", $TAC/$tot);
        $TTA = sprintf("%.17f", $TTA/$tot); $TTT = sprintf("%.17f", $TTT/$tot); $TTG = sprintf("%.17f", $TTG/$tot); $TTC = sprintf("%.17f", $TTC/$tot);
        $TGA = sprintf("%.17f", $TGA/$tot); $TGT = sprintf("%.17f", $TGT/$tot); $TGG = sprintf("%.17f", $TGG/$tot); $TGC = sprintf("%.17f", $TGC/$tot);
        $TCA = sprintf("%.17f", $TCA/$tot); $TCT = sprintf("%.17f", $TCT/$tot); $TCG = sprintf("%.17f", $TCG/$tot); $TCC = sprintf("%.17f", $TCC/$tot);

        $GAA = sprintf("%.17f", $GAA/$tot); $GAT = sprintf("%.17f", $GAT/$tot); $GAG = sprintf("%.17f", $GAG/$tot); $GAC = sprintf("%.17f", $GAC/$tot);
        $GTA = sprintf("%.17f", $GTA/$tot); $GTT = sprintf("%.17f", $GTT/$tot); $GTG = sprintf("%.17f", $GTG/$tot); $GTC = sprintf("%.17f", $GTC/$tot);
        $GGA = sprintf("%.17f", $GGA/$tot); $GGT = sprintf("%.17f", $GGT/$tot); $GGG = sprintf("%.17f", $GGG/$tot); $GGC = sprintf("%.17f", $GGC/$tot);
        $GCA = sprintf("%.17f", $GCA/$tot); $GCT = sprintf("%.17f", $GCT/$tot); $GCG = sprintf("%.17f", $GCG/$tot); $GCC = sprintf("%.17f", $GCC/$tot);

        $CAA = sprintf("%.17f", $CAA/$tot); $CAT = sprintf("%.17f", $CAT/$tot); $CAG = sprintf("%.17f", $CAG/$tot); $CAC = sprintf("%.17f", $CAC/$tot);
        $CTA = sprintf("%.17f", $CTA/$tot); $CTT = sprintf("%.17f", $CTT/$tot); $CTG = sprintf("%.17f", $CTG/$tot); $CTC = sprintf("%.17f", $CTC/$tot);
        $CGA = sprintf("%.17f", $CGA/$tot); $CGT = sprintf("%.17f", $CGT/$tot); $CGG = sprintf("%.17f", $CGG/$tot); $CGC = sprintf("%.17f", $CGC/$tot);
        $CCA = sprintf("%.17f", $CCA/$tot); $CCT = sprintf("%.17f", $CCT/$tot); $CCG = sprintf("%.17f", $CCG/$tot); $CCC = sprintf("%.17f", $CCC/$tot);


return "$AAA $AAT $AAG $AAC $ATA $ATT $ATG $ATC $AGA $AGT $AGG $AGC $ACA $ACT $ACG $ACC $TAA $TAT $TAG $TAC $TTA $TTT $TTG $TTC $TGA $TGT $TGG $TGC $TCA $TCT $TCG $TCC $GAA $GAT $GAG $GAC $GTA $GTT $GTG $GTC $GGA $GGT $GGG $GGC $GCA $GCT $GCG $GCC $CAA $CAT $CAG $CAC $CTA $CTT $CTG $CTC $CGA $CGT $CGG $CGC $CCA $CCT $CCG $CCC";
}
################################
sub categoriza_rnasize {
      my ($tamanho_seq) = @_;
      my $conteudo_processado;
      if ($tamanho_seq <= '100') {
#       $conteudo_processado = "105:1 106:0 107:0 108:0"; #SEM \n PORQUE É FORNECIDO DIRETO Á LINHA FINAL SVM
       $conteudo_processado = "85:1"; #SEM \n PORQUE É FORNECIDO DIRETO Á LINHA FINAL SVM
      }
      elsif (($tamanho_seq > '100') && ($tamanho_seq <= '400')) {
#       $conteudo_processado = "105:0 106:1 107:0 108:0";
       $conteudo_processado = "86:1";
      }
      elsif (($tamanho_seq > '400') && ($tamanho_seq <= '900')) {
#       $conteudo_processado = "105:0 106:0 107:1 108:0";
       $conteudo_processado = "87:1";
      }
      elsif ($tamanho_seq > '900') {
#       $conteudo_processado = "105:0 106:0 107:0 108:1";
       $conteudo_processado = "88:1";
      }

 return($conteudo_processado);
}
##########################

