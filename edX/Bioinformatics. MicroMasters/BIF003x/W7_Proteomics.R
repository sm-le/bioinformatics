# Week 7 Proteomic
# Sungmin Lee

# Objectives
# Describes the concept of Microarray analysis
# Demonstrates how to perform a T test on our adult vs fetal brain
# Discusses the significance of False Discovery Rate (FDR) correction
# Discusses the text output file and how to plot at different expression levels

#######################
# Microarray Analysis #
#######################
# one of the most common microarray analysis is expression analyses
# where the spot on the chip represent different genes, and mRNA isolated from an experimental sample
# which are treated by 1) the RNA binds to its complimentary DNA, 2) probes on the RNA that can be used to detect binding levels

library(affy)
library(limma)

# Load all CEL files
affy.data <- ReadAffy()

# The dataset was created using the mas5 probe, so we need to normalize the data to compensate difference between runs
eset.mas5 <- mas5(affy.data)

# Retrieve normalized matrix and store it as an expression set
exprSet.nologs <- exprs(eset.mas5)

# rename column names
colnames(exprSet.nologs) <- c("brain.1","brain.2","fetal.brain.1","fetal.brain.2","fetal.liver.1","fetal.liver.2","liver.1","liver.2")

# log transform to represent fold change in expression
exprSet <- log(exprSet.nologs,2)

# export exprSet
write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")

# Not all microarray spots will be populated with data, so
# it is useful to "Absent/Present" check.
data.mas5calls <- mas5calls(affy.data) # generates a vector with the A/P values for the dataset
data.mas5calls.calls <- exprs(data.mas5calls) # converts into an expression matrix containing A or P for each tissue and gene combination. 

# limma package to read images of brain microarray
brain.fetalbrain.2color <- read.maimages("brain.fetalbrain.2color.dat.txt",columns=list(G="brain.1",R="fetal.brain.1", Gb="bg1", Rd="bg2"))

# normalize color # Loess normalization uses a form of local-weighted regression to standardize the data
brain.fetalbrain.2color.loess <- normalizeWithinArrays(brain.fetalbrain.2color, method='loess')

# two panel graph result
par(mfrow=c(1,2)) # tells R that our plot window should have 1 row and 2 columns
plotMA(brain.fetalbrain.2color) # plot original image data
plotMA(brain.fetalbrain.2color.loess) # corrected data

# calculate the log base 2 ratio of expression for comparison of binary expression level between conditions.
# first get mean of each dataset
brain.mean <- apply(exprSet[,c("brain.1","brain.2")],1,mean)
fetal.brain.mean <- apply(exprSet[,c("fetal.brain.1","fetal.brain.2")],1,mean)
liver.mean <- apply(exprSet[,c("liver.1","liver.2")],1,mean)
fetal.liver.mean <- apply(exprSet[,c("fetal.liver.1","fetal.liver.2")],1,mean)
# (1, means) is the MARGIN for the command, which indicates that it is the rows within each column that the mean is calculated for

# log transformation of the expression values
brain.fetal.to.adult <- fetal.brain.mean - brain.mean
liver.fetal.to.adult <- fetal.liver.mean - liver.mean

# make dataframe and export
all.data <- cbind(exprSet, brain.mean, fetal.brain.mean, liver.mean, fetal.liver.mean, brain.fetal.to.adult, liver.fetal.to.adult)
write.table(all.data, file = "Microarray_ALL.txt", quote=F, sep='\t')

# perform t-test on adult vs fetal brain
dataset.1 <- exprSet[1,c("brain.1","brain.2")]
dataset.2 <- exprSet[1,c("fetal.brain.1","fetal.brain.2")]

t.test.gene.1 <- t.test(dataset.1, dataset.2, "two.sided")

# apply a t-test on relevant column
brain.p.value.all.genes <- apply(exprSet, 1, function(x){t.test(x[1:2],x[3:4]) $p.value})
# taking exprSet with a margin of 1 (to look at paired item in a row) and applying t-test on merges of column 1 and 2 and compares them to columns 3 and 4

# on liver data set
liver.p.value.all.genes <- apply(exprSet, 1, function(x) { t.test(x[5:6], x[7:8]) $p.value } )

# A/P test and filter out uninformative genes.
AP <- apply(data.mas5calls.calls, 1, paste, collapse="")
length(AP)
# applies the paste command across the data.mas5calls.calls dataset at the row level (margin =1) and places the A or P values into a new variable, AP

# any gene with all As are uninformative
genes.present = names(AP[AP != "AAAAAAAA"])
length(genes.present)

# make subset
exprSet.present <- exprSet[genes.present,] # only getting rows with data

# Look at FDR
brain.raw.pvals.present <- brain.p.value.all.genes[genes.present]
liver.raw.pvals.present <- liver.p.value.all.genes[genes.present]

brain.fdr.pvals.present <- p.adjust(brain.raw.pvals.present, method="fdr")
liver.fdr.pvals.present <- p.adjust(liver.raw.pvals.present, method="fdr")

# sort gene by p-value
brain.fdr.pvals.present.sorted <- brain.fdr.pvals.present[order(brain.fdr.pvals.present)]
liver.fdr.pvals.present.sorted <- liver.fdr.pvals.present[order(liver.fdr.pvals.present)]

# look at first 10
brain.fdr.pvals.present.sorted[1:10]
liver.fdr.pvals.present.sorted[1:10]

# output with very high p-value
brain.DE.probesets <- names(brain.raw.pvals.present[brain.raw.pvals.present < 0.01])
liver.DE.probesets <- names(liver.raw.pvals.present[liver.raw.pvals.present < 0.01])

# Put log 2 ratio
brain.DE.log2.ratios <- all.data[brain.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]
liver.DE.log2.ratios <- all.data[liver.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]

# export
write.table(brain.DE.log2.ratios, "brain.DE.log2.ratios.txt", sep="\t", quote=F)
write.table(liver.DE.log2.ratios, "liver.DE.log2.ratios.txt", sep="\t", quote=F)

# Note
# How to test different gene expression
# Parametric test
# 1. two sample t-test
# 2. Welch t-test
# Non-parametric tets
# 1. Permutation t-test
# 2. Wilcoxon rank sum test

# Permutation t-test
# - resampling method, if the labels are exchageable under the null hypothesis, then the resulting tests yield exact significance levels
# - arrange the combined set of data in ascending order, receiving a rank equal to average
# - U statistic, U=n1*n2+{n1*(n1+1)/2} -T (n1, n2 size of first and second sample)

# When testing so many gene, it is easy to fasely call many genes under significance
# define False Positive Gene is important in this case.
# Benjamin Hochberg Method, method for choosing how many genes to call significantly differentially expressed
# - p <= x*a/m
# 1. p -> largest p-value
# 2. a -> FDR rate e.g) 0.05
# 3. x -> number of genes to call significant
# 4. m -> total number of genes

# Clustering
# to find inner data structure according to similarity, group similar observation
# define the distance by dissimilarity, function of data
# 1. hierarchical clustering
# 2. centroid based clustering (partition)

#################
# Visualization #
#################

# a plot of expression levels is equally valuable in terms of analytics

# single plot mode
par(mfrow=c(1,1))

# create plot
x.data <- all.data[, "brain.mean"]
y.data <- all.data[, "fetal.brain.mean"]

plot(x.data, y.data, main = "Log2 expression in fetal brain (n=2) vs adult brain (n=2)", xlab="brain", ylab="fetal brain", col="blue", cex=0.5)
abline(0,1)

# ratio intensity plot, commonly used for genomic plot or bland-altman
# display the difference between measurements between two samples by transforming that data onto a log ratio (M) scale and average (A) scale.
brain = all.data[,"brain.mean"]
fetal.brain = all.data[,"fetal.brain.mean"]
A = (fetal.brain + brain)/2
M = (fetal.brain-brain)

plot(A, M, main="MA plot of fetal brain vs brain", pch=19, cex=0.2, col="red")

# Volcano plot 
# plot significance versus fold-change 
expression.plus.pvals=cbind(exprSet.present, brain.raw.pvals.present, brain.fdr.pvals.present, liver.raw.pvals.present, liver.fdr.pvals.present)

log2.ratios = expression.plus.pvals[,"brain.1"] - expression.plus.pvals[,"fetal.brain.1"]
p.values = expression.plus.pvals[,"brain.raw.pvals.present"]

plot(log2.ratios, -log(p.values, 10) )

# generate a side-by-side comparative plot of the volcano plot we just created and one where we instead don't do a log base 10 transform
par(mfrow=c(1,2))
plot(log2.ratios, p.values)
# To read a volcano plot, you look for the points that are furthest away from the bottom center – those points represent the genes that are most deferentially expressed.

# Summary
# The concept of Microarray analysis
# How to perform a T test on our adult vs fetal brain data
# The significance of False Discovery Rate(FDR) correction
# How to generate a text output file and to plot at different expression levels







