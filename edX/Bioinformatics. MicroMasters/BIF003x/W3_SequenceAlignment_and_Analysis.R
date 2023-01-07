# Week 3
# Sequence Alignment & Analysis

# Objectives
# Demonstrates protein sequence alignment for the Histone H1 protein taken from 5 different species
# Illustrates how to create dotplots from two sequences called dotPlots
# Describes the different types of multiple sequence alignment algorithms and parameters
# Demonstrates how to export multiple sequence alignments
# Highlights the different features and methods of phylogenetic reconstruction

# Sungmin Lee
# Pairwise Sequence Alignment

# Introduction
# One of the most basic functions in sequence analysis is a sequence alignment.
# A sequence alignment takes two or more sequences that are suspected of being similar and "lines them up," one on top of the other, matching similar residues (bases) and inserting gaps (dashes) where necessary.
# Sequence alignments can also include "mismatches" – residues or bases that do not match but are orthologous – meaning that they are known to be in the same position in both sequences despite not being identical.
# Sequence alignments are generally the first step in any comparative analysis.

# Fetch bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library(Biostrings)
library(seqinr)

# Load DNA sequences
prokaryotes <- read.fasta(file="prok.fasta",seqtype="DNA")

# split files into individual sequences for pairwise alignment
seq1 <- as.character(prokaryotes[[1]])
seq1 = paste(seq1,collapse="")
seq2 <- as.character(prokaryotes[[2]])
seq2 = paste(seq2,collapse="")

# align seq 1 and seq2
pairalign <-  pairwiseAlignment(pattern=seq2,subject=seq1)

# export alignment out of R as a FASTA file
pairalignString = BStringSet(c(toString(subject(pairalign)), toString(pattern(pairalign))))

writeXStringSet(pairalignString,"aligned.txt",format="FASTA")

# sequence comparison visuals with Dot Plot
coxgenes <- read.fasta(file="cox1multi.fasta",seqtype="AA") # aa file is imported as DNA dot plot can be messy
cox1 <- as.character(coxgenes[[1]])
cox2 <- as.character(coxgenes[[2]])

dotPlot(cox1,cox2,main="Human vs Mouse Cox1 Dotplot")
# To improve signal to noise ratio, consider
# 1. wsize <- moving window of size wsize
# 2. wstep <- number of step between windows
# 3. nmatch <- control how many residue/bases within the window must match

# reduce noise
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# only first 100 aa
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Dot Plot interpretation
# off diagonal line show inverted transposition.
# to interpret this, place x and y coordinate and start writing.

#In a local alignment, you are matching your query with a substring (fragment) of your subject sequence. 
#In a global alignment you perform an end to end alignment between the two.
#You may end up with a lot of gaps in global alignment if the sizes of query and subject are dissimilar).
#Local alignments may also have gaps.

# Multiple Sequence Alignment

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("msa")

library(msa)

# Use readAAStringSet and readDNAStringSet from Biostrings
coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

#StringSets are collections of processed sequences and can be of several types:

#DNAStringSet
#RNAStringSet
#AAStringSet
#BStringSet

#The first three are typed to only allow a specific type of sequence, e.g. DNA, RNA, or AA.
#A BStringSet can contain any type of sequence.
#Using a typed StringSet enables other analyses to "know" what type of sequence to expect – meaning we don't have to declare the type in many calls.

# Mutltiple Sequence Alignment with StringSet as an input
coxAligned <- msa(coxAA)

prokAligned <- msa(prokDNA)

# display alignment

print(prokAligned,show="complete")

# Options in MSA
msa(prokDNA, "ClustalW")
msa(prokDNA, "ClustalOmega")
msa(prokDNA, "Muscle")

# Parameters in MSA
# 1. cluster <- which clustering method
# 2. gapOpening <- determines gap opening penalty
# 3. gapExtension <- determines gap extension penalty
# 4. Maxiters <- maximum number of iteration
# 5. SubstitutionMatrix <- specify a substitution matrix 
# 6. Order type <- aligned vs input, dna or rna or protein
# 7. Verbose <- progress message or not

# Two ways of exporting MSA, fasta or phylip

# Exporting as FASTA
prokAlignStr = as(prokAligned,"DNAStringSet") # has to match input alignment type
writeXStringSet(prokAlignStr,file="prokAligned.fasta")

coxAlignStr = as(coxAligned,"AAStringSet")
writeXStringSet(coxAlignStr,file="coxAligned.fasta")

# Exporting as Phylip
write.phylip(coxAligned,"coxAligned.phylip") # do not have to convert the alignment to a stringset

# Most simply put, MSAs are used to look at the similarities – and differences – between groups of sequences.

# Phylogenetic Reconstruction

# 1. distance
# 2. maximum parsimony
# 3. maximum likelihood

# prokaryotic alignment to seqinr format
prokAligned2 <- msaConvert(prokAligned, type="seqinr::alignment")

# distance matrix using seqinr and ape
prokdist <- dist.alignment(prokAligned2,"identity")
prokdist

# Neighbour-joining method for distance phylogenetic tree
library(ape)
prokTree <- nj(prokdist)
plot(prokTree)

# Maximum parsimony with phangorn
library(phangorn)
prokAligned3 <- msaConvert(prokAligned, type="phangorn::phyDat")
ParsTree <- pratchet(prokAligned3)
plot(ParsTree)

# Maximum likelihood
# Calculate likelihood with pml function
fit <- pml(prokTree, prokAligned3)
# improve tree with JC model
fitJC <- optim.pml(fit,model="JC",rearrangement="stochastic") # can use other model like K80
plot(fitJC)

# Bootstrapping
# method of subsampling data
bootstrapped <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control=pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bootstrapped, p=50, type="p")

# Summary
#Protein sequence alignment for the Histone H1 protein taken from 5 different species
#How to create dotplots from two sequences called dotPlots
#The different types of multiple sequence alignment algorithms and parameters
#How to export multiple sequence alignments
#The different features and methods of phylogenetic reconstruction


