#!/bin/sh

bin_dir=$1
blast_ver="ncbi-blast-2.7.1+-x64-linux.tar.gz"
blast_path="ncbi-blast-2.7.1+/bin/"
cpc2="CPC2-beta.tar.gz"
cpc2_path="CPC2-beta/bin/CPC2.py"
d2p="dna2pep-1.1.tar.gz"
d2p_path="dna2pep-1.1/dna2pep.py"
libsvm="libsvm-3.18.tar.gz"

cd $bin_dir

######################NEEDS SUDO PERMISSIONS
## update & upgrade #
#sudo apt-get update
#sudo apt-get upgrade

echo "INSTALLING ANNOCRIPT PRE-REQUISITES"
#MySQL
echo "INSTALLING MySQL Server and Client"
sudo apt-get install -y mysql-server
apt-get install mysql-client
#PERL
echo "INSTALLING PERL"
sudo apt-get install -y perl
#BioPERL
echo "INSTALLING BIOPERL"
sudo apt-get install -y bioperl
#PERL Libraries
echo "INSTALLING Additional PERL libraries"
sudo apt-get install -y libterm-readkey-perl
sudo apt-get install -y libparallel-forkmanager-perl
#PYTHON
echo "INSTALLING PYTHON"
sudo apt-get install -y python
#R
echo "INSTALLING R"
sudo apt-get install -y r-base
######################NEEDS SUDO PERMISSIONS


##CPC2
echo "INSTALLING CPC2"
gzip -dc $cpc2 | tar xf -
cd CPC2-beta
export CPC_HOME="$PWD"
cd libs/libsvm
gzip -dc $libsvm | tar xf -
cd libsvm-3.18
make clean && make


#Return to DL_PROGRAMS
cd ../../../../
#BLAST
echo "INSTALLING BLAST "
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/$blast_ver

cp $blast_ver $bin_dir
cp $cpc2 $bin_dir
cp $d2p $bin_dir
cd $bin_dir
gzip -dc $cpc2 | tar xf -
gzip -dc $d2p | tar xf -
gzip -dc $blast_ver | tar xf -
rm $blast_ver


echo "##############################PROGRAM LINKS #################################"
echo "Please use the following links for the programs into the configuration file:"
echo blastPath = $bin_dir/$blast_path
echo lncRNAPredProgPath = $bin_dir/$cpc2_path
echo dna2pepPath = $bin_dir/$d2p_path
echo "##############################PROGRAM LINKS #################################"


exit 1
