#Creates a BAR plot for the enrichments


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


##Installs temporarily missing packages
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {
#	install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
  tmp.install.packages("RColorBrewer")
}
#else{print("RColorBrewer already installed")}


library(RColorBrewer)


args <- commandArgs(trailingOnly = TRUE)

#Stores in the matrix the table give in input
table <- read.table(args[1], header=TRUE, sep="\t", quote="", stringsAsFactors = FALSE, dec=".")

#Creates these two lists for the values
values <- head(table$occurrences,20)
descs <- head(table$description,20)


#creates an image to save
jpeg(args[2], width = 1024, height = 768, pointsize = 12, quality= 100, unit = "px") 

#Creates a color vector to use to color the bars
color_vector <- colorRampPalette(brewer.pal(9,"Set1"),bias=1 )( 20 )

#Set the margins of the plot to leave a space on the right
par(mar = c(2, 2, 2, 50), xpd = NA)

#Creates the plot
barplot(values, col=color_vector, main=args[3], xlab=args[4], ylab=args[5])


#Creates the legend
#'x' defines the position respect to the plot and 'inset' define the distance from the margins
legend(legend=descs, bty="n", fill = color_vector, cex=0.9, , x= "right",inset=-2)

#Closes the image
dev.off()



