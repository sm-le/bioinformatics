# Week 8 R for bioinformatics
# Sungmin Lee
par(mfrow=c(1,1))
# Objectives
# Discusses a high-throughput method of sequencing RNA – essentially transcriptome sequencing
# Demonstrates how to use Bioconductor to examine RNAseq data
# Highlights the different ways to analyze quick statistics and how to manipulate different data sets
# Demonstrates how to generate and manipulate different plots under quality control

##########################
# RNAseq Data Processing #
##########################

# RNAseq is high-throughput method of sequencing RNA - essentially transcriptome sequencing.
# rapidly replacing microarray assays when it comes to expression studies

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)

# read RNAseq count
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt",stringsAsFactors = FALSE) # contain gene count info
sampleinfo <- read.delim("Sampleinfo.txt") # contain the metadata about our samples

# filter data to contain only counts
countdata <- seqdata[,-(1:2)]

# pushing first column of seqdata to the rownames of countdata
rownames(countdata) <- seqdata[,1]

# truncate column names
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)

# filter lowly expressed genes
myCPM <- cpm(countdata)
thresh <- myCPM > 0.5

# check filter result
table(rowSums(thresh))

# filter out lowest value
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]

# check if all has correct CPM number
plot(myCPM[,1],countdata[,1])
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)

##########################
# RNAseq Quality Control #
##########################

# convert our truncated data from counts.keep to the dgelist format required for the EdgeR library
y <- DGEList(counts.keep)

# check if there are any discrepancies between the samples
barplot(y$samples$lib.size, names=colnames(y),las=2) # the las=2 command rotates axis text.
title("Barplot of library sizes") 

# check further
logcounts <- cpm(y, log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# the data is roughly evenly distributed amongst the samples

# Another good GC plot is MDSplot, stands for Multidimensional Scaling, a visual representation of a principle components analysis
par(mfrow=c(1,2))
levels(sampleinfo$CellType)
col.cell <- c("red","blue")[sampleinfo$CellType]

# plot MDS
plotMDS(y,col=col.cell)
legend("topleft",fill=c("blue","red"),legend=levels(sampleinfo$CellType))
title("Cell type")

# color plot by condition
levels(sampleinfo$Status)
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"), legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# import corrected data
sampleinfo <- read.delim("SampleInfo_Corrected.txt")
par(mfrow=c(1,2))
col.cell <- c("blue","red")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("blue","red"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"), legend = levels(sampleinfo$Status), cex=0.8)
title("Status")

################################
# Analysis - MDS Dimesionality #
################################

# We can also us pch= to change the character plotted
par(mfrow=c(1,1))
char.celltype <- c("X", "O") [sampleinfo$CellType]
plotMDS(y,col=col.status,pch=char.celltype,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))
title("MDS Plot for first two dimension")

#  add the dim command and plot alternate dimensions
plotMDS(y,dim=c(3,4), col=col.status,pch=char.celltype,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))
title("MDS Plot for dimensions 3 and 4") 

# Hierarchical clustering
# large datasets the best way to visualize the results of clustering are heatmaps
# gplots contains a command, heatmap.2, that will both perform a hierarchical clustering and generate a heatmap

# Compute variance in each row
var_genes <- apply(logcounts, 1, var)

# Limit to 500 most variant genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

# subset our logcounts by choosing only the genes that show up in select_var
highly_variable_lcpm <- logcounts[select_var,]

# summary of the contents of highly_variable_lcpm
dim(highly_variable_lcpm) # 500 genes, and 12 values for expression

# color scheme
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette) # generate a function rather than discrete color

# perform clustering and heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples", ColSideColors = col.cell,scale="row")
# using 50 shades of color sorted in reverse order from our morecols function
# The trace="none" command instructs heatmap.2 not to draw connecting lines through columns or rows

# save img
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()
# dev.off() terminate output redirect
# or we can do dev.copy2pdf() without switching to redirect
# advantage of higher resolution 

# Summary
# The high-throughput method of sequencing RNA - essentially transcriptome sequencing
# How to use Bioconductor to examine RNAseq data
# The different ways to analyze quick statistics and how to manipulate different data sets
# How to generate and manipulate different plots under quality control