#tables_2_barplot.R
#Copyright (C) <2014>  <Francesco Musacchia>
#Realized for inclusion in Annocript 0.2

#
#Creates a BAR plot using a table given in input which has two columns:
# 'occurrences' and 'descriptions' 

#Gets the parameters in input
args <- commandArgs(trailingOnly = TRUE)

#Table in input with abundances
in_table <- args[1]

#name of the image to print
img_name <- args[2]

#Name of the figure
main_name <- args[3]

#X axis name
x_name <- args[4]

#Y axys name
y_name <- args[5]

#Maximum number of elements to show in the plot
top_to_show <- as.numeric(args[6])

#Define color palette
col.vec<-c("#A6CEE3","#87BAD8","#69A7CD","#4B94C3","#2C80B8","#3084AE","#519BA5","#72B29C","#93C992","#AFDD88","#92CF72","#76C15D","#59B348","#3DA533","#4F9F3B","#7C9D54","#A99C6C","#D69B84","#FA9493","#F47877","#EF5B5B","#E93E3F","#E42123","#E73429","#ED593C","#F27F4E","#F8A461","#FDBB68","#FDAC4F","#FE9E36","#FE8F1D","#FE8104","#F58827","#E99357","#DD9F87","#D1AAB7","#C2A8D1","#AC8EC3","#9773B6","#8159A8","#6B3F9B","#886499","#A99099","#CBBB99","#ECE799","#F7EE8D","#E5C874","#D4A35A","#C27E41","#B15928")

###################################COMMENTED
#This code is commented because I decided to use a defined color palette instead that
#ColorBrewer. This last should be installed before but sometimes R does not permit.
#You can uncomment this code if you want and instead comment the line
#with the definition of color_vector

##creates a temporary directory and install there the package
##It wil be erased at the next system restart
#tmp.install.packages <- function(pack, dependencies=TRUE) {
#  path <- tempdir()
#  ## Add 'path' to .libPaths, and be sure that it is not
#  ## at the first position, otherwise any other package during
#  ## this session would be installed into 'path'
#  firstpath <- .libPaths()[1]
#  .libPaths(c(firstpath, path))
#  install.packages(pack, dependencies=dependencies,repos="http://cran.us.r-project.org", lib=path)
#}


###Installs temporarily missing packages
#if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {
##	install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
#  tmp.install.packages("RColorBrewer")
#}
##else{print("RColorBrewer already installed")}

#library(RColorBrewer)
#color_vector <- colorRampPalette(brewer.pal(9,"Set1"),bias=1 )( as.numeric(args[7]) )
#Creates a color vector to use to color the bars

#COMMENT THIS IF YOU WANT TO USE COLORBREWER
color_vector <- col.vec[sample(50,top_to_show)]


#Stores in the matrix the table give in input
table <- read.table(in_table, header=TRUE, sep="\t", quote="", stringsAsFactors = FALSE, dec=".")

#Creates these two lists for the values
values <- head(table$occurrences,top_to_show)
descs <- head(table$description,top_to_show)

#creates an image to save
jpeg(img_name, width = 1024, height = 768, pointsize = 12, quality = 100, unit = "px") 


#Set the margins of the plot to leave a space on the right
par(mar = c(2, 2, 2, 50), xpd = NA)

#Creates the plot
barplot(values, col = color_vector, main = main_name, xlab = x_name, ylab = y_name)


#Creates the legend
#'x' defines the position respect to the plot and 'inset' define the distance from the margins
legend(legend = descs, bty = "n", fill = color_vector, cex = 0.9, , x = "right", inset = -2)

#Closes the image
dev.off()



