# R Test Script
# Sungmin Lee
# Nov 11th, 2019
# Test script to learn R

# To learn and utilize variables and data in R
count <- 0

# To make complex data structure, a vector

primes <- c(1,3,5,7,11)

# c is telling R to create a list

Names <- c("Bob", "Ted", "Carol", "Alice")
Truth <- c(TRUE, FALSE)

# DataFrame is a table with row and column

organism <- c("Human","Chimpanzee","Yeast")
chromosomes <- c(23,24,16)
multicellular <- c(TRUE, TRUE, FALSE)

OrganismTable <- data.frame(organism, chromosomes, multicellular)

# Data Structure Simple Command

OrganismTable$organism
OrganismTable$organism[2]

# Write data.frame to table in csv format
write.table(OrganismTable, file = "MyData.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

# Using count to report how many organisms have more than 20 chromosomes

# direct counting
if (OrganismTable$chromosomes[1] > 20) count=count+1
if (OrganismTable$chromosomes[1] > 20) count=count+1
if (OrganismTable$chromosomes[1] > 20) count=count+1

# For loop counting

count<-0
for (val in OrganismTable$chromosomes) {
  if (val > 20) 
    count = count+1
}
print(count)

# plotting data - organism table

barplot(OrganismTable$chromosomes)

# import ggplot2 for better aesthetic

library(ggplot2)
library(dplyr)


rawdata <- read.csv("Week_1_Plotdata.csv",header=TRUE)

# plotting rawdata

rawdata %>%
  ggplot(aes(x=Subject,y=a))+
  geom_point()

# melt our data into different form

library(reshape2)

# reshape to cut and paste desired data from rawdata

melted = melt(rawdata, id.vars="Subject", measure.vars=c("a","c","d","e","f","g","j","k"))

myPlot <- melted %>%
  ggplot(aes(x=variable,y=value,col=Subject,group=Subject)) +
  geom_point() +
  geom_line() +
  xlab("Sample") +
  ylab("# Observed") +
  ggtitle("Some observation I made in the lab")

myPlot

ggsave(filename="ParcuPrevalance.pdf",plot=myPlot)
