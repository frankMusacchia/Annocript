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
--------------------------
PORTRAIT - Prediction of transcriptome ncRNA by Ab initio Methods
--------------------------
Version for peer reviewing purposes

Authors: -

Webserver version may be found at:
http://bioinformatics.cenargen.embrapa.br/portrait

--------------------------
Installation instructions
--------------------------
Download software package, then on a Linux terminal go to the path where package is located.
Type:
gunzip -c portrait-1.1.tar.gz | tar xvf -
 
This should create a new directory, PORTRAIT. Enter this directory:
cd PORTRAIT
Now install LIBSVM by first unpacking it:
gunzip -c libsvm-2.84.tar.gz | tar xvf -
Then enter LIBSVM directory, and compile it:
cd libsvm-2.84
make
Afterwards, PORTRAIT is ready to be installed.
Go to PORTRAIT directory and type:
./portrait-1.1.pl
You will be asked to inform where ANGLE and CAST are located.
These are included in the package so typically you will find them in PORTRAIT directory/cast-linux/ and /angle/.

All programs are automatically tested by PORTRAIT on-the-fly during installation.
If you manage to reach the end of installation process without any error messages or warnings, PORTRAIT will be ready to use.

Usage:
./portrait-1.1.pl -i <transcripts file> -s <libsvm> -c <cast> -a <angle> [options]
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
        -h : This help.

--------------------------
Known Bugs and Issues
--------------------------
Q1. If you get an error message like this one:
"ANGLE-linux64DP: cannot execute binary file"
A1. Most likely your library libstdc++ is either outdated or missing. Confirm by executing
either ANGLE-linux32DP or ANGLE-linux64DP alone. ANGLE requires libstdc++ version 6. 
PORTRAIT cannot work properly until this problem is fixed.

Q2. Some of my input sequences are missing in the output file.
A2. This was a bug present on PORTRAIT early versions: it ignored some sequences if user FASTA input
had a small number of sequences (about 5 sequences or less). This is mainly because PORTRAIT was from
beginning crafted for working with large-scale data. This bug was fixed after version 1.0,
but should you encounter this problem, just add some dummy sequences to the end of your input.

