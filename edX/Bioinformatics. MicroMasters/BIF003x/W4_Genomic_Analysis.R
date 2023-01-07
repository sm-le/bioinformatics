# Genomic Analysis
# Sungmin Lee
# Objective

# 1. Demonstrate how to sort all of our outputs numerically using the sort function in Biostrings
# 2. Highlight the different ways features of common to genes like ORF, promoters, CPG islands, and splice sites
# 3. Discuss how to identify START and STOP sites • Describes the different ways of modifying the code to plot results
# 4 .Illustrate how to comparisons between genomes using Comparative Genomics

# Genomic Analysis
# Predictive: predict where genes and other sequence elements may lie on genome
# Comparative: find sequence regions that are homologous to known sequence region
# Both approaches are prone to errors

# Predictive anaylsis -> gene prediction
# look for common features: 1. ORF, 2. Promoters, 3. CpG islands, and 4. Splice Sites

# ORF Finding - example of string search
# 1. Length of the ORF: minimum length must be established
# 2. Reading frames: there are 6 potential reading frames for any given sequence
# Start codon: ATG, Stops are TAG, TAA,TGA
# We can use biostrings or built-in ORF finder.

library(Biostrings)
library(seqinr)

AB003468 <- readDNAStringSet("AB003468.fasta")
AB003468 <- as.character(AB003468)

# use matchPattern from biostrings to search for any pattern
# find start codon
matchPattern("ATG",AB003468)

sequence <- AB003468

# make search more portable, create start and stop codon variables
start_codon <- "ATG"
stop_codons <- c("TGA","TAA","TAG")

# we make 4 result vars, create vars to store results in
start_pos <- c()
revstart_pos <- c()
stop_pos <- c()
revstop_pos <- c()

# Forward matches
matches <- matchPattern(start_codon,sequence)
start_pos <- c(start_pos, start(matches))

# Reverse matches
revmatches <- matchPattern(reverseComplement(DNAString(start_codon)),sequence)
revstart_pos <- c(revstart_pos,start(revmatches))

# Sort Results
start_pos <- sort(start_pos)
revstart_pos <- sort(revstart_pos, decreasing = TRUE)

# Finding stop codons
for (codon in stop_codons){
  matches <- matchPattern(codon,sequence)
  stop_pos <- c(stop_pos, start(matches))
  revmatches <- matchPattern(reverseComplement(DNAString(codon)),sequence)
  revstop_pos <- c(revstop_pos, start(revmatches))
}

stop_pos <- sort(stop_pos)
revstop_pos <- sort(revstop_pos, decreasing = TRUE)

# set ORF threshold for minimum size
k <- 150 # set minimum ORF to 50 aa.
lengths <- vector(mode="numeric")

# Forward frame
stop_pointers <- c(0,0,0) # hold the locations of the STOPS in each reading frame
count <- 0 # number of ORF we found

# Find ORF above threshold in forward frame

for (current_start in start_pos){ # loop through each possible starting position
  frame <- (current_start%%3) + 1 # keep track of the frame
  stop_pointer <- stop_pointers[frame] # most recent stop
  
  if (stop_pointer <= length(stop_pos)  # This loop tell the current stop pointer has to be smaller than the length of stop position
      && (stop_pointer == 0 # if the stop pointer is still zero
          || stop_pos[stop_pointer] < current_start)) { # or if stop_pos of the current position is less than the current start position
    
    stop_pointer <- stop_pointer+1 # we increment stop_pointer and continue
    
    while ((stop_pointer <= length(stop_pos)) # while stoppointer is less than stop_pos length
           && ((stop_pos[stop_pointer] <= current_start)  # current stop pos is less than current start
               || (((stop_pos[stop_pointer]%%3) + 1) != frame)) # stop_pos in the current frame # basically scanning from start to end until stop
          ){
      stop_pointer <- stop_pointer+1
    }
    stop_pointers[frame] <- stop_pointer
    
    if (stop_pointer <= length(stop_pos)) {
      if ((stop_pos[stop_pointer]+2-current_start +1) > k) { # check length of our potential ORF, by stop_pos[stop_pointer]
        count <- count+1 # print out all current variables
        print(count)
        print("Frame:")
        print(frame)
        print("Start:")
        print(current_start)
        print("Stop:")
        print(stop_pos[stop_pointer])
        print("Length:")
        lengths <- c(lengths, (stop_pos[stop_pointer]+2-current_start+1)) # we are adding length of ORF to vector legnths
        print(stop_pos[stop_pointer]+2-current_start+1)
        print("Sequence:")
        print(subseq(sequence, current_start, stop_pos[stop_pointer]+2))
      }
    }
  }
}


#  Find ORF above threshold in reverse frame

revstop_pointers <- c(0,0,0)

for (current_revstart in revstart_pos) {
  current_revstart <- current_revstart + 2
  frame <- (current_revstart%%3) + 1
  revstop_pointer <- revstop_pointers[frame]
  if (revstop_pointer <= length(revstop_pos) && (revstop_pointer == 0
                                                 || revstop_pos[revstop_pointer] > current_revstart)) {
    
    revstop_pointer <- revstop_pointer +1
    
    while ((revstop_pointer <= length(revstop_pos))
           && ((revstop_pos[revstop_pointer] + 2 >= current_revstart)
               || ((((revstop_pos[revstop_pointer] + 2)%%3) + 1) != frame))) {
      revstop_pointer <- revstop_pointer + 1
    }
    revstop_pointers[frame] <- revstop_pointer
    
    if (revstop_pointer <= length(revstop_pos)) {
      if ((current_revstart - revstop_pos[revstop_pointer])+1 > k) {
        count <- count + 1
        print(count)
        print("Frame:")
        print(-frame)
        print("Start:")
        print(current_revstart)
        print("Stop:")
        print(revstop_pos[revstop_pointer] + 2)
        print("Length:")
        lengths <- c(lengths, (current_revstart - revstop_pos[revstop_pointer]))
        print(current_revstart - revstop_pos[revstop_pointer])
        print("Sequence:")
        print(subseq(sequence,revstop_pos[revstop_pointer],current_revstart))
      }
    }
  }
}

lengths <- sort(lengths)
lengths

# Plotting the result

barplot(lengths) # barplot
plot(density(lengths)) # density plot

bins <- seq(0,1000,50) # kernel density plot
hist(lengths, breaks=bins, col="red", xlim = c(0,1000))

# Comparative genomics

# Using Biomart get ensembl database
library(biomaRt)
ensembl = useMart("ensembl")

# To access available dataset list
listDatasets(ensembl)

# Load Chimpanzee and Gorilla dataset
Chimpanzee <- useDataset("ptroglodytes_gene_ensembl", mart=ensembl)
Gorilla <- useDataset("ggorilla_gene_ensembl", mart=ensembl)
Human <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Pull list of available attributes
ChimpanzeeAttributes <- listAttributes(Chimpanzee)
GorillaAttributes <- listAttributes(Gorilla)
HumanAttributes <- listAttributes(Human)


# Retrieve unique gene identifier
ChimpanzeeGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Chimpanzee) # getBM retrieve data from ensembl server
GorillaGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Gorilla)
HumanGenes <- getBM(attributes = c("ensembl_gene_id","gene_biotype"),mart=Human)

# Parse Gene name and type into different vectors
ChimpanzeeGeneNames <- ChimpanzeeGenes[[1]]
ChimpanzeeGeneTypes <- ChimpanzeeGenes[[2]]
GorillaGeneNames <- GorillaGenes[[1]]
GorillaGeneTypes <- GorillaGenes[[2]]

# Summarize type using table
ChimpanzeeTypeTable <- table(ChimpanzeeGeneTypes)
GorillaTypeTable <- table(GorillaGeneTypes)

# Compare number of protein coding in gorilla vs chimpanzee
ChimpanzeeProteins <- ChimpanzeeTypeTable["protein_coding"]
GorillaProteins <- GorillaTypeTable["protein_coding"]

# Number of gorilla orthologs in Chimpanzee
ChimpGGNum <- getBM(attributes = "ggorilla_homolog_ensembl_gene", mart=Chimpanzee)
ChimpanzeeGG <- getBM(attributes = c("ensembl_gene_id", "ggorilla_homolog_ensembl_gene"), mart=Chimpanzee)

# Number of chimpanzee orthologs in Gorilla
GorillaCNum <- getBM(attributes = "ptroglodytes_homolog_ensembl_gene", mart = Gorilla)
GorillaC <- getBM(attributes = c("ensembl_gene_id","ptroglodytes_homolog_ensembl_gene"), mart = Gorilla)

# Summary

# How to sort all of our outputs numerically using the sort function in Biostrings
# The different features common to genes like ORF, promoters, cpG islands and splice sites
# How to identify START and STOP sites
# The different ways of modifying the code to plot results
# How to compare between genomes using Comparative Genomics