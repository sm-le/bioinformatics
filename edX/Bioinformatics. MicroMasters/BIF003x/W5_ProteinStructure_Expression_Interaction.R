# Week 5 Protein Structure, Expression, and Interaction
# Sungmin Lee

library(gplots)

# Objective
# 1. Demonstrates how to install gplots and the various features of distribution plots
# 2. Highlights the different ways of analyzing data with multiple conditions using Linear Discriminate Analysis (LDA)
# 3. Discusses how to install BioConductor packages and describes the functions of litG library
# 4. Describes the methods of visualizing graphs and and extracting subgraphs

################################
# Protein Analysis - Mass Spec #
################################

# MALDI TOF - general background
# Larger fragments go shorter distances; smaller fragments go further.
# By looking at all possible fragment sizes generated, software can assemble a list of all probable proteins in a sample of mixed proteins – as well as their quantity.

# Essentially two properties: Protein ID and Protein quantities

# Distribution plot of peptides from a MALDI-TOF
peptides.txt <- read.table("peptidefrags.txt",header=FALSE)
peptides <- as.vector(peptides.txt$V1)
hist(peptides,breaks=400)

# Another QC test between three different versions of mass-spec analysis
mascot.txt <- read.table("mascot.txt",header=FALSE)
xtandem.txt <- read.table("xtandem.txt",header=FALSE)
protpro.txt <- read.table("protpro.txt",header=FALSE)

mascot <- as.vector(mascot.txt$V1)
xtandem <- as.vector(xtandem.txt$V1)
protpro <- as.vector(protpro.txt$V1)

# Quick note about list
# 1. Lists can mix data types.
# 2. Lists do not need to be defined as a class type.
# 3. Most importantly for this exercise – lists don't need every element (in this case, our vectors) to be the same size.

# Combine MS result 
combinedMSdata <- list(Mascot=mascot, XTandem=xtandem, ProtPro=protpro)

# Use venn diagram for summary result
venn(combinedMSdata)

################
# LDA Analysis #
################

library(timeSeries)
library(MASS)
library(rgl)
library(ggplot2)

# Load ms data
Dataset <- read.csv("ms.csv",header=TRUE, na.strings="NA",dec=".",strip.white=TRUE)

# Limit dataset to actual protein data
RawData <- Dataset[,2:14]

# To fill empty columns
filledcols = colSds(RawData) != 0.0 # wherever the value is not zero, creates a vector called filledcols listing those columns
RawData <- RawData[,filledcols] # replacing with only those columns in filledcols

# specify LDA class
test1.lda <- lda(Dataset$X1 ~ . , data=Dataset) # creates LDA probabilities for each class in X1
test1.lda.values <- predict(test1.lda,Dataset) # we can generate actual LDA values by testing against the data
# You will always have N-1 LDA scores for a dataset which contains N conditions, e.g. if you had 6 conditions you would end up with 5 LDA score columns as the output.

# visualize LDA result
x <- test1.lda.values$x[,1]
y <- test1.lda.values$x[,2]
# The LDA values are always listed in order of largest contribution, so while you may have more columns to choose from, plotting the first two is always a good move. 

# distinguish which point from which condition
class <- Dataset$X1
plotdata<-data.frame(class,x,y)

centroids <- aggregate(cbind(x,y)~class, plotdata, mean) # Set center of cluster

CentroidDistances <- dist(centroids, method = "euclidean", diag=TRUE, upper=FALSE, p=2) # generate a distance matrix from the dataframe centroids using straight-line
attr(CentroidDistances, "Labels") <- centroids$class

# plot
plot1 <- plotdata %>% ggplot(aes(x,y,color=factor(class))) +
  geom_point(size=3) +
  geom_point(data=centroids,size=7)+
  geom_text(data=centroids, size=7, label=centroids$class, color="black") +
  ggtitle("LDA of Conditions 1-3") +
  geom_text(aes(label=Dataset$X1),hjust=0,vjust=0,color="black")
plot1

# save LDA centroid distance data
write.csv(as.matrix(CentroidDistances, file = "centroiddistances.csv"))

#################################
# Protein - Protein Interaction #
#################################

library(graph)
library(Rgraphviz)
library(RBGL)
library(yeastExpData) # https://www.nature.com/articles/ng776z nat-gen article about the dataset


# From the yeastExpData library let's load the litG dataset
data("litG") # load data

# Extract all nodes from litG
litGnodes <- nodes(litG) # names of all proteins

adj(litG, "YFL039C") # which proteins are connected to each other

# invidual connected components
connectedComponents <- connectedComp(litG)

# Pull out subgraph3 
component3 <- connectedComponents[[3]]
subgraph3 <- subGraph(component3, litG)

# Plot
subgraph3plot <- layoutGraph(subgraph3, layoutType="neato")
renderGraph(subgraph3plot)

# Number of degrees for each protein
numdegrees <- graph::degree(litG)
graph::degree
numdegrees <- sort(numdegrees)

meandeg <- mean(numdegrees)

hist(numdegrees, col="red", main=paste("Degree Distribution - Protein Interaction in litG with a mean of " , meandeg))

###########
# Summary #
###########
# 1. How to install gplots and the various features of distribution plots
# 2. The different ways of analyzing data with multiple conditions using Linear Discriminate Analysis (LDA)
# 3. How to install BioConductor packages and the functions of litG library
# 4. The methods of visualizing graphs and extracting subgraphs


####################################
# Access Mass Spec raw data from R #
####################################

library(BiocInstaller)
biocLite("MSnbase")
library(msdata) # example ms data 
library(mzR) # fundamental infrastructure for all word dealing with mass spec
library(MSnbase)

fls <- proteomics(full.names=TRUE)
basename(fls)
fl <- fls[2]
fl

rw <- openMSfile(fl)
rw # contains 565 spectra, not loaded but can be accessed

sp1 <- spectra(rw, 1) 
spl <- spectra(rw)

hd <- header(rw)
head(hd)

mse <- readMSData(fl, mode = "onDisk")
mse[[1]]
fData(mse)
