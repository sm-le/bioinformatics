# R script for Basic DNA analysis
# Sungmin Lee

# Working with DNA sequence

library(seqinr)

# read fasta file

cox1 <- read.fasta(file="cox1.fasta",seqtype="AA")

length(cox1)
seq1 <- cox1[1]

# retrieve DNA sequence from Genbank

library(ape)

AB003468 <- read.GenBank("AB003468",as.character="TRUE")

# Save GenBank sequence to fasta file

write.dna(AB003468,file="AB003468.fasta",format="fasta",append=FALSE, nbcol=6,colsep=" ",colw=10) 

# 1. export can be done by seqinr or ape libraries, 2. should be save and written in plain text format

# retrieve sequence using NCBI Entrez
library("rentrez")

entrez_search(db="nucleotide",term="human superoxide dismutase")

# Sequence Statistics

# Convert sequence into a simple string
CloningVector <- AB003468[[1]]

count <- count(CloningVector,1) # Count number of single nucleotide, it is word count summary
# count(CloningVector,2) for dinucletoide and count(CloningVector,3) for trinucleotide

# To compute GC contect
GC <- GC(CloningVector) # GC affect melting temperature of DNA

# It would be interesting to see if GC value changes over a length of nucleotide - GC content by window
GCwindow <- seq(1,length(CloningVector)-200,by=200) # -200 to avoid going beyond the sequence limit

n <- length(GCwindow) # number of chunk of dna to analyze

Chunks <- numeric(n) # blank vector with the same number of blank space

# For loop to compute GC chunk 
for (i in 1:n){
  chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)] # +199 for 1-200, 201-400 operation
  chunkGC <- GC(chunk)
  print(chunkGC)
  Chunks[i] <- chunkGC # to add GC content into a vector
}

plot(GCwindow, Chunks, type="b", xlab="Nucleotide start position", ylab="GC content")

# Custom functions in R, a way to define a piece of code that you can invoke. This won't be built in and need to be called everytime
# Write Custom function in R that will calculate GC contect given only input sequence and window size

slidingwindowGCplot <- function(windowsize, inputseq){
  GCwindow <- seq(1, length(inputseq)-windowsize, by=windowsize)
  # Find length of GCwindow
  n <- length(GCwindow)
  # Make a blank vector which length equilivalent to n
  Chunks <- numeric(n)
  
  for (i in 1:n){
    chunk <- inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    Chunks[i] <- chunkGC
  }
  plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position", ylab = "GC content", main=paste("GC Plot with windowsize ", windowsize)) # past command let us combine multiple items into the same title
}

# Test custom function
slidingwindowGCplot(100,CloningVector)

# Protein Sequence Statistics

library(Peptides)

# amino acid composition of the sequence
aaComp(cox1[1])

# Alipatic index of protein sequence
aIndex(cox1) # indicator of thermostability

# predict net charge of a protein
charge(cox1)

# specify sequence in the command for net charge
charge(seq="FLPVLAG",pH=7, pKscale="EMBOSS")

# calculate hydrophobicity
hydrophobicity(cox1[1])

'
Week 2 summary
- Demonstrated how to read and write FASTA
- Highlighted the different simple sequence statistics and retrieving sequences from NCBI
- Discussed the features of GC content determination and plotting
- Described the functions of word content determination and plotting
'
