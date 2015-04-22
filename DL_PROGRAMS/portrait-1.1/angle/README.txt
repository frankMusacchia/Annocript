****************************************************************************************
                       ANGLE: Anti Noise Gene Locator Engine
                  Copyright (C) Kana Shimizu, All rights reserved, 2007
****************************************************************************************

1. USAGE

Copy the folder into your PC and execute ANGLE from a console as follows.

ANGLE [OPTION]

OPTIONS:
        -i: Input filename (default name is sample.txt.)
        -d: Paramater directory (default directory is ./param-human)
        -w: Display width (default width is 60 bases.)
        -h: Output frame layout to a html file.
        -s: Output score of each base to a file.

----
ex) ANGLE -w 70 -h out.html -s score.cvs -d paramater_directory -i inputFasta.txt

2. Input File Style
ANGLE accepts fasta format WITH UPPER-CASE LETTERS as input files. If your file includes lower-case lettres, please uppercase them.

3. Output Style

ANGLE has three types of output.

3.1	Default
Input sequences are displayed with coding sequence information.
- Coding sequence are shown between "{"  "}".
- Detected insertion errors are displayed between "(" ")"
- Detected deletion errors are displayed as "X".

3.2	HTML option
- Coding sequence is shown in red color.
- Detected insertion errors are striked-through with blue color.
- Detected deletion errors are displayed as "X" with blue color.
HTML option also shows graphical output as follows.

CGGCAGCTCGGATGCGGCGGACGGCGGGACA 
-------------------------------  <- first frame (begins from the first base)
-------------------------------  <- second frame (begins from the second base)
-----------********************  <- third frame (begins from the third base)

"-" shows the codon in the frame is not included in coding sequence. "*" shows the codon in the frame is included in coding sequence.
For example, "ATG" which is located in third frame is included in coding sequence, but next three bases "TGC" which is located in first frame is not included in coding sequence.

3.3	Score option
ANGLE assigns each base a score, which shows how much the base is likely to be included in coding sequence. If the score for a base is larger than 0, the base is predicted to be included in coding sequence. This score file is written in CVS format, so that the text information can be easily converted to other spreadsheet softwares like EXCEL.

4. Copyright notice
The program is free for scientific purposes, as long as its content is not modified. Please contact the author, if you are planning to use the software for commercial purposes. The software must not be further distributed without prior permission of the author. 

5. REFERENCE

Shimizu, K., Adachi, J., and Muraoka, Y. (2006) ANGLE: a sequencing errors resistant program for predicting protein coding regions in unfinished cDNA. J Bioinform Comput Biol 4, 649-664

