#Annocript 0.2 - Output

##How *Annocript* manages folders and files


When you run *Annocript* the first time two folders are created in the working directory you created 
(e.g. ann_works) (Here we will write folders in capital letters):

- JOBS: contains all the analysis data you will run and get
- DATA: contains all the file downloaded for the database and the index files for blast

Inside the JOBS folder you will see all the sessions you saved as folders. Suppose that you session is S1. The folder S1 will have the following content:
- DATA: a folder containing a set of files that *Annocript* creates during the execution.        
  - outHashFile.hash: an hash that will contain all the output from the programs to be used successively for fast execution
  - configHash: the hash containing all the parameters of the programs
-OUTPUT: contains all the output by *Annocript* and the programs. Thus it contains:
  - BLASTX: a folder containing results from blastx (blastp);
  - BLASTN: a folder containing results from blastn (tblastn);
  - RPSTBLASTN: a folder containing results from rpstblastn (rpsblast);
  - portraitOut: output from Portrait execution;
  - dna2pepOut: output from the dna2pep execution;
  - S1_uniref_2014_04_ann_out.txt: is the complete tabular output of *Annocript*;
  - S1_uniref_2014_04_filt_ann_out.txt: is the reduced tabular output of *Annocript* (without the sequence);
  - S1_uniref_2014_04_lncRNAseqs.fasta: is a fasta file with all the putative lncRNAs identified by *Annocript*. See below to understand;
  - S1_uniref_2014_04_orf_info.fasta: is a fasta file with all the longest ORF as extracted by dna2pep. It is a fasta file containing in the header;
    - the name of the transcript
    - the strand and the frame as found from dna2pep and blastx (blastp)
    - the longest ORF
  - S1_uniref_2014_04_codingSeqs.fasta: a fasta file with all the coding sequences as found by *Annocript*. See below to understand.ì

- GFF: this is a folder containing the output from blast programs given in GFF3 format: rpstblastn_out.gff, blastn_out.gff, portait_out.gff, blastxUniref_out.gff, blastxSprot_out.gff, dna2pep_out.gff. They are overwritten if the GFF are created again.
- STATS: This folder contains:
a set of files with percentages of presence of GO classes, pathways and organisms/taxonomies: GO_bp_percentages, GO_mf_percentages, GO_cc_percentages, pwl1_percentages, pwl2_percentages, pwl3_percentages, closer_os_table.txt; ATGC percentages: ATGC_Percentages.txt; The percentage is computed as the sum of all the occurrences of the term divided per the number of annotated sequences;
a HTML file (comprehensive of a folder with images) with all the statistics and plots coming out from the files of percentages;
- S1_UNIREF_STATSPAGE: is the folder coming with the HTML file containing all the needed images;
       
- LOG: it is a folder containing:
  - a log file for *Annocript*. It comes with date and time to be specific of a particular run; 
  - a file with times occurred to run *Annocript* that has date and time both files have date and time.
  - R log, coming out from the generation of plots in R.
- S1_SLICE: this folder contains you transcriptome divided in equally distributed bases. These  files are created to be ran parallely with rpsblast.
- portrait_S1.fasta: it is the fasta file created ad-hoc to use in Portrait. Since some stuffs in the nucleotidic sequences are not accepted by the software, it should be cleaned;
- S1.fasta: the transcriptome inserted by the user;
- config_user.txt: the file with the configuration used.


##Log files

Apart from the main log file about the run of *Annocript* and the file with times, a log files is created by *Annocript* during the construction of DB.
uniprotkb_headParsNAValues.log: this log file is written while the UniProt_kb/UniRef table construction is running. This table will contain a set of data fetched from the headers of the corresponding databases from UniProt. If a tag does not contain anything, *Annocript* will write in this log the corresponding row and tag which does not appear. You can consider to look at this file to see what is the progress in the reading of the UniProt database;
R.log: will be in the LOG folder and will contain the output of the R script for the plots.

##What does *Annocript* gives in the tabular output

###Columns of the output file with Uniprot_KB 
(UniProt Knowledgebase containing both results from SwissProt and TrEMBL)

*Annocript* gives in output a text file with tab separated fields. They are the following:
	
```
                TranscriptName: is the name of the transcript as given in your FASTA file
                TransLength: is the length in nucleotides (or aminoacids) of the transcript
                Count: this is the number of reads that have been used from the assembler program to build this transcript
                HSPNameSP: this is the first HSP result (lowest e-value) as given from the BLASTx output against Swiss-Prot.
                HSPLengthSP: is the corresponding length in nucleotides of the HSP result
                HSPEvalueSP: is the corresponding e-value assigned to the HSP
                HITLengthSP: is the length of the first HIT as given from blastx output
                QCoverageSP: this is how much of the query (transcript) is covered from the HSP (See below)
                HCoverageSP: this is how much the HSP covers the HIT (See below)
                DescriptionSP: description of the HSP
                EnzymeIds: enzyme identifier corresponding to the HSP. (See below)
                EnzymeDescs: enzyme descriptions of the  EnzymeIds (See below)
                HSPNameTR:  this is the first HSP result (lowest e-value) as given from the blastx output against Trembl
                HSPLengthTR: is the corresponding length 
                HSPEvalueTR: is the corresponding e-value assigned to the HSP
                HITLengthTR: is the length of the HIT as given from blastx output
                QCoverageTR: this is how much of the query (transcript) is covered from the HSP (See below)
                HCoverageTR: this is how much the HSP covers the HIT (See below)
                DescriptionTR: description of the HSP
                OSName: organism corresponding to the result among SwissProt and TrEMBL with lowest e-value 
                BPId: Biological processes ids corresponding to the lowest e-value result between SwissProt and TrEMBL (separated by ]---[) (See below)
                BPDesc: Biological processes descriptions with same sorting as in BPId (See below)
                MFId: molecular functions ids with same sorting as in BPId (See below)
                MFDesc: molecular descriptions  with same sorting as in BPId (See below)
                CCId: cellular components ids  with same sorting as in BPId (See below)
                CCDesc: cellular components descriptions  with same sorting as in BPId (See below)
                CDName: all domains names coming out from the output of rpsblast
                CDStartEnd: domains start-end coming out from the rpsblast
                CDEvalue: domains e-values 
                CDDesc: domains descriptions
                domBPId: Biological processes ids corresponding to the results of rpsblast against pfam domains (each pfam may contain more GO terms separated with ';'. Pfam results are in turn separated by ]---[) (See below)
                domBPDesc: Biological processes descriptions with same sorting as in domBPId (See below)
                domMFId: molecular functions ids with same sorting as in domBPId (See below)
                domMFDesc: molecular descriptions with same sorting as in domBPId (See below)
                domCCId: cellular components ids with same sorting as in domBPId (See below)
                domCCDesc: cellular components descriptions with same sorting as in domBPId (See below)
                OtherNCName: name of the ribosomial or other ncRNA sequences coming from the HSP with lowest e-value output of blastn (See below)
                OtherNCEvalue: e-value of the ribosomial or other ncRNA sequences coming from the output of blastn (See below)
                OtherNCDesc: descriptions of the ribosomial or other ncRNA sequences coming from the output of blastn (See below)
                LongOrfLength: length of the longest ORF as found from dna2pep
                LongOrfStrand: strand of the longest ORF as found from dna2pep
                LongOrfFrame: frame of the longest ORF as found from dna2pep
                ProbToBeNonCoding: probability that the sequence is non-coding coming from Portrait
                lncRNA4Annocript: Annocript's heuristic about the sequence. If it is 1, it is a lncRNA. See below.
                Sequence: the aminoacidic (peptidic) sequence of the transcript
```

The filtered version of the output does not contain the sequences.

###Columns of output file the with UniRef 

It contains both results from SwissProt and UniRef. The table is the same but now the GO terms are always associated with Uniref result. And the following names change to "Uf" for UniRef.

                HSPNameUf: this is the HSP with lowest e-value as given from the blastx output against UniRef
                HSPLengthUf: is the corresponding length of the HSP
                HSPEvalueUf: is the corresponding e-value assigned to the HSP
                HITLengthUf: is the length of the HIT as given from blastx output
                QCoverageUf: this is how much of the query (transcript) is covered from the HSP (See below)
                HCoverageUf: this is how much the HSP covers the HIT (See below)
                DescriptionUf: description of the HSP 
                Taxonomy: taxonomy corresponding to the result of UniRef
                
               
The filtered version of the output does not contain the sequences.

**lncRNA4Annocript**

This field is 1 when *Annocript* thinks that the sequence is a lncRNA. You will see the lncRNA sequences inside a FASTA file with all of them.

The sequence is non coding if:
- there is no homology result with Proteins, Domains and other non-coding RNAs;
- the longest orf length is less or equal than NCORFLength (default: 100);
- the length of the query is equal or more than NCSeqLength (default: 200);
- ProbToBeNonCoding is more or equal to NCThresh (default: 0.95).
               
Where NCORFLength, NCSeqLength, NCThresh are modifiable parameters in the user configuration file.


**Coding sequences**

The file with the coding sequences contains sequences which have at least a homology result among SwissProt, TrEMBL/UniRef or domains or they should have the longest ORF bigger than NCORFLength (default: 100).

**QCoverage and HCoverage**

They are computed as follows:
- QueryCoverage = (length of query participating in alignment minus gaps / query length)*100;
- HitCoverage = (length of hit participating in alignment minus gaps / hit length)*100;

**EnzymeIds**

This value is extracted by using the file at the following link:
    
    enzymeDBLink = ftp://ftp.expasy.org/databases/enzyme/enzyme.dat

The database built by *Annocript* contains a table where each enzyke identifier is associated with the corresponding SwissProt id. In the output table we take the one associated to "HSPNameSP".

**GO classes**

User can choose if wants to display GO terms associated to proteins or domains. Biological process (BP), molecular functions (MF) and cellular components (CC) come out from the GO ids associated to the UniRef/UniProt_kb ids. We extracted the mapping for SwissProt, UniRef and TrEMBL proteins from the file at:      

    GODBLink = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

While mapping for Pfam domains are taken at:
    
    http://geneontology.org/external2go/pfam2go

**Other NC**

A search against other noncoding RNAs (rRNA, tRNA, snRNA, snoRNA, miRNA) is performed with a blastn search against an integrated database of Rfam and NCBI Refseq ribosomal RNAs.


----------------------------------------

If you get some error during the installation or the running of Annocript please see the FAQ page of the user guide or ask on the google group: https://groups.google.com/forum/#!forum/annocript
