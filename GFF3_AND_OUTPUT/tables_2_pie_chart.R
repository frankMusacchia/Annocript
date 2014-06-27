#Creates a pie chart for the enrichments


###################################USEFUL

#creates a temporary directory and install there the package
#It wil be erased at the next system restart
tmp.install.packages <- function(pack, dependencies=TRUE) {
  path <- tempdir()
  ## Add 'path' to .libPaths, and be sure that it is not
  ## at the first position, otherwise any other package during
  ## this session would be installed into 'path'
  firstpath <- .libPaths()[1]
  .libPaths(c(firstpath, path))
  install.packages(pack, dependencies=dependencies,repos="http://cran.us.r-project.org", lib=path)
}

##Installing missing packages
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {
#	install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
  tmp.install.packages("RColorBrewer")
}
#else{print("RColorBrewer already installed")}

library(RColorBrewer);
#Create a PIE graphic for the occurrences

args <- commandArgs(trailingOnly = TRUE)

#Store in a matrix the table given in input
table <- read.table(args[1],header=TRUE, stringsAsFactors=FALSE,sep="\t",quote='')

#take only the first 20 values from the occurrences column
values <- head(table$occurrences,20)

#take only the first 20 values from the description column
descs <- head(table$description,20)

jpeg(args[2], width = 1024, height = 768, pointsize = 12, quality= 100, unit = "px") 



par(mar=c(2,2,2,40),xpd=TRUE,bty="n")
color_vector <- colorRampPalette(brewer.pal(9,"Set1"),bias=1 )( 20 )

#Create the pie
pie(values, labels = '',col=color_vector,radius=-1, main=args[3])

box()

legend(1.1,1.0,legend=descs, bty="n", fill = color_vector);

#hist(values, labels = descs, main=args[3])
dev.off()
