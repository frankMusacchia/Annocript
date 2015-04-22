#-------------------------------------------------
# PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------
#R CMD BATCH --no-save --no-restore '--args <job to do> <annocript_filt_out>  <out folder> <name>' ann_plot_stats.R
# - Takes in input the annocript output
# - a parameter that indicates what to produce, 
# - the folder where to put the output
# - a name to give to the output file

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)
#The samples consiedered are 4 - please comment whenever you will not use all
length(args)
if (length(args)<3){
	stop("USAGE: R CMD BATCH --no-save --no-restore '--args <job to do> <annocript_filt_out> <out folder> <name>'  ann_plot_stats.R\n
	1: transcripts lengths histogram\n2: distribution of lengths barplot\n3: distribution of coverage barplot ")
}

#Folder where all is taken and output goes
job = args[1]
job

#Annocript output path
ann = args[2]
ann

#output folder path
out_folder = args[3]
out_folder

#Name for the output
name = args[4]
name




ann.out = read.delim(file=ann,sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)

# 1: transcripts lengths histogram
if ((job == 1) || (job == 'all')){
	jpeg(file=paste(out_folder,paste(name,'transc_lengths.jpeg',sep='_'),sep="/"), width = 680, height = 480, units = "px", pointsize = 12,
     bg = "white")
	#freq = (max(ann.out$TransLength)-min(ann.out$TransLength))/100
	freq = 100
	hist(ann.out$TransLength,breaks = freq,  main="Histogram of Lengths", xlab="Transcripts lengths",col="purple")
	#axis(1, at = seq(min(ann.out$TransLength), max(ann.out$TransLength), by = freq*4), las=2)
	dev.off()
}

# 2: distribution of lengths barplot
if ((job == 2) || (job == 'all')){
	jpeg(file=paste(out_folder,paste(name,'distrib_lengths.jpeg',sep='_'),sep="/"), width = 680, height = 480, units = "px", pointsize = 12,
     bg = "white")
	
	boxplot(ann.out$TransLength, main="Distribution of Lengths")
	dev.off()
}

# 3: distribution of query coverage barplot
if ((job == 3) || (job == 'all')){
	jpeg(file=paste(out_folder,paste(name,'distrib_q_coverages.jpeg',sep='_'),sep="/"), width = 680, height = 480, units = "px", pointsize = 12,
     bg = "white")
	
	if ("QCoverageUf" %in% colnames(ann.out)){
		boxplot(as.numeric(ann.out$QCoverageSP),as.numeric(ann.out$QCoverageUf),main="Distribution of Query Coverages",names=c("SwissProt","UniRef"),
		col=c("blue","purple"),ylab = "Percentage")
	}
	if ("QCoverageTR" %in% colnames(ann.out)){
		boxplot(as.numeric(ann.out$QCoverageSP),as.numeric(ann.out$QCoverageTR),main="Distribution of Query Coverages",names=c("SwissProt","TrEMBL"),
		col=c("blue","purple"),ylab = "Percentage")
	}
	
	dev.off()
}

# 4: distribution of hit coverages barplot
if ((job == 4) || (job == 'all')){
	jpeg(file=paste(out_folder,paste(name,'distrib_h_coverages.jpeg',sep='_'),sep="/"), width = 680, height = 480, units = "px", pointsize = 12,
     bg = "white")
	
	if ("HCoverageUf" %in% colnames(ann.out)){
		boxplot(as.numeric(ann.out$HCoverageSP),as.numeric(ann.out$HCoverageUf),main="Distribution of Hit Coverages",names=c("SwissProt","UniRef"),
		col=c("blue","purple"),ylab = "Percentage")
	}
	if ("HCoverageTR" %in% colnames(ann.out)){
		boxplot(as.numeric(ann.out$HCoverageSP),as.numeric(ann.out$HCoverageTR),main="Distribution of Hit Coverages",names=c("SwissProt","TrEMBL"),
		col=c("blue","purple"),ylab = "Percentage")
	}

	dev.off()
}

# 5: longest ORF lengths histogram
if ((job == 5) || (job == 'all') && is.numeric(ann.out$LongOrfLength)){
	jpeg(file=paste(out_folder,paste(name,'lORF_lengths.jpeg',sep='_'),sep="/"), width = 680, height = 480, units = "px", pointsize = 12,
     bg = "white")
	#freq = (max(ann.out$LongOrfLength)-min(ann.out$LongOrfLength))/100
	freq = 100
	hist(ann.out$LongOrfLength,breaks = freq,  main="Histogram of Longest ORF Lengths", xlab="lORFs lengths",col="blue")
	#axis(1, at = seq(min(ann.out$TransLength), max(ann.out$TransLength), by = freq*4), las=2)
	dev.off()
}
