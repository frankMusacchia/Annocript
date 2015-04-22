#!/usr/bin/perl -w
use warnings;
use strict;

# RTARRIAL MAR/2007

###############
# Afold feeder
###############
# O programa de predicao de estrutura secundaria de RNA, Afold, aceita
# apenas uma sequencia por rodada. Esse script, entao, lê um arquivo multifasta
# e a cada nova sequencia ele fornece ao Afold, gerando um arquivo multioutput.

# Supondo que o multifasta de entrada já tenha sido filtrado de seqs > 28.000 letras,
# mas se tiver alguma sequencia assim, ela será eliminada antes de ser fornecida ao Afold

#O script só funciona se for rodado dentro do mesmo diretório onde está o Afold!

# ESPECIFIQUE O DIRETORIO do PROGRAMA AFOLD::
my $afold_dir = '/home/rtarrial/Afold_RNA';

my $usage = "<arquivo multifasta> <arquivo de saida multioutput>";
if ($#ARGV < 0) { die "$0 $usage";}
my $infasta = $ARGV[0];
my $output = $ARGV[1];

my $temp_in = "$afold_dir/afold_feeder.temp";

existe_arquivo($infasta);
nao_existe_arquivo($output);
nao_existe_arquivo($temp_in);
nao_existe_arquivo("$afold_dir/Afold.res");

open(OUT, ">$output") or die "\nErro ao salvar para $output: $!\n";
my $notfirsttime = '0';
my $sequencia;
my $header;

open(FASTA, $infasta) or die "\n$!";

while (my $linha = <FASTA>) {
  chomp $linha;
  if ($linha =~ /#/) {next;} #ignorando partes comentadas
  elsif ($linha =~ /^\s*$/) { next;} #ignorando espaços em branco

  elsif (($linha =~ /^\>/) && ($notfirsttime)) { #Se está lendo uma nova seq, primeiro processe a seq anterior..
    my $tamanho = length($sequencia);
    if ($tamanho > 28000) {
     print "\nSeq: $header foi excluida - possui $tamanho letras";
    }
    else {
           open(TEMP, ">$temp_in") or die "\nErro ao salvar para temp $temp_in: $!\n";
           chomp $header;
           print OUT "$header\n";
           print TEMP "$header\n$sequencia"; #criando arquivo temp de entrada para Afold
           close TEMP;
           `$afold_dir/Afold $temp_in`;
           unlink $temp_in; #apagar arquivo temp
           open(RES, "$afold_dir/Afold.res") or die "\nErro no arquivo salvo pelo Afold: $! \n";
           while (<RES>) {
            print OUT $_;
           }
          # print OUT "---\n"; #Separar duas predicoes
           unlink("$afold_dir/Afold.res"); #apagar output do Afold
    }
  $sequencia = '';
  $header = $linha;
  }
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

    my $tamanho = length($sequencia);
    if ($tamanho > 28000) {
     print "\nSeq: $header foi excluida - possui $tamanho letras";
    }
    else { open(TEMP, ">$temp_in") or die "\nErro ao salvar para temp $temp_in: $!\n";
           chomp $header;
           print OUT "$header\n";
           print TEMP "$header\n$sequencia"; #criando arquivo temp de entrada para Afold
           close TEMP;
           `$afold_dir/Afold $temp_in`;
           unlink $temp_in; #apagar arquivo temp
           open(RES, "$afold_dir/Afold.res") or die "\nErro no arquivo salvo pelo Afold: $! \n";
           while (<RES>) {
            print OUT $_;
           }
     #      print OUT "---\n"; #Separar duas predicoes
           unlink("$afold_dir/Afold.res"); #apagar output do Afold
           }

close FASTA;
print "\nMultiOut salvo com sucesso em $output.\n";

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


