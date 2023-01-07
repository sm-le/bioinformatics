# Sungmin Lee
# Algorithm for DNA sequencing
# Week 1: DNA sequencing, strings, and matching

###################
# Why study this? #
###################
'''
- DNA sequencing has become inexpensive, and there are numerous data
- There are many applications that studies various areas.
1. Understanding key algorithm is key to understand where they succeed and failed.
2. is the first step to figure out how to do better.
'''
###################################
# DNA sequencing past and present #
###################################
'''
1. 1st Generation sequencing, require manual intervention of human.
- Used to be a dominant sequencing technology.

2. 2nd Generation sequencing, or next generation sequencing
- Massive parallel sequencing
'''
###########################################
# Genomes as strings, reads as substrings #
###########################################
'''
Genome is book of recipe in analogy, which is writeen in A, C, G, and T.

If we read DNA from up to down, we can write DNA molecule as a string.

DNA sequencers are not good at reading a long stretch of DNA. Instead they read short sequences.
They read randomly selected substrings or called sequencing reads (100nt).
'''
##########################################
# String definitions and Python examples #
##########################################
'''
String S is a finite sequence of characters
Characters are drawn from alphabet sigma, where sigma = {A,C,G,T}
|S| = number of characters in S

s = 'ACGT'
-> len(s)
s[2]
-> 'G'
s='AACC'
t='GGTT'
s+t
-> 'AACCGGTT'
s='AACCGGTT'
s[2:6]
-> 'CCGG'
s[0:6]
-> 'AACCGG'
'''
#################
# Strings basic #
#################

# define seq
seq = 'ACGT' # or double quote is fine
# a character in index 1 of string
seq[1]
# to get a length of character
len(seq)
# empty string
e = ''
# concatenate string
seq1='CCAA'
seq2='GGTT'
print(seq1+seq2)
# concatenate with join function
seq = ['A','C','G','T']
print(','.join(seq))
# random draw
import random
random.choice('ACGT') # if want same result random.seed(n)
# generate sequence from random sample
seq = ''
for _ in range(10): # don't save number in variable so used underscore (_)
    seq += random.choice('ACGT')
print(seq)

seq = ''.join([random.choice('ACGT') for _ in range(10)]) # other form

# get substrings
seq[0:3] # or seq[:3]

# or get suffix
seq[7:]

############################
# Manipulating DNA strings #
############################

# find longest common prefix
def longestCommonPrefix(s1, s2):
    i = 0
    while i < len(s1) and len(s2) and s1[i] == s2[i]:
        i+=1
    return s1[:i]

# find if sequences matche
def match(s1, s2):
    if not len(s1) == len(s2):
        return False
    for i in range(len(s1)):
        if not s1[i] == s2[i]:
            return False
    return True

s1 == s2 # simply can check if they match

# get reverse compliment

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

####################################
# Downloading and parsing a genome #
####################################
#!wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa
#curl https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa

# read genome from fasta file
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# count frequency of each base
counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for base in genome:
    counts[base] += 1

# easier count
import collections
collections.Counter(genome)

#######################
# How DNA gets copied #
#######################

'''
complementary strands are separated and act as templates
DNA polymerase integrate individual nucleotide to the strand
'''

#########################################
# How second generation sequencer works #
#########################################

'''
1. cut input DNA into snippets
2. deposit snippets on a slide, there are billions of template scattered across the plate
3. add DNA polymerase, and bases
4. each base has terminator, polymerase can only add one complimentary base
5. take photo, photo signal emitted from each bases
6. remove terminators and repeat
'''

######################################
# Sequencing errors and base quality #
######################################

'''
cluster forms in next generation sequencing for signal amplification. In this cluster, un-terminated base can be added and add additional bases (out-of-sync)
As a result, multiple light can be emitted from one cluster

base color will evaluate the quality of the sequencing and report base quality (in a probability term).
Q = -10 * log(10)p
Q = base quality
p = probability that base call is incorrect

how base quality conveys uncertainty?
Estimate p for example, non-orange light/total light (simple assumption)
'''

####################################
# Sequencing reads in FASTQ format #
####################################

'''
In a FASTQ format
first line: Name of the read (experience, sequencing machine, where on the slide)
second line: Sequence of bases
thrid line: ignore
fourth line: sequence of base quality

Each base quality is ASCII-encoded version of Q=-10log(10)p
Usual ASCII encoding is Phred+33
take Q, rounded to integer, add 33, convert to character
'''

##############################
# Work with sequencing reads #
##############################

# read FASTQ file
def readFASTQ(filename):
    sequences = []
    qualities = []
    with open(filename,"r") as f:
        while True:
            f.readline().rstrip()
            seq = f.readline().rstrip()
            f.readline().rstrip()
            qual = f.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

seqs, quals = readFASTQ('SRR835775_1.first1000.fastq')

# what quality score is most common?
# First convert Phred+33 to quality score

def phred33ToQ(qual):
    return ord(qual) - 33

# Create histogram from quality score
def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

h=createHist(quals)
import matplotlib.pyplot as plt
%matplotlib inline
plt.bar(range(len(h)),h)
plt.show()

###############################
# Analyzing reads by position #
###############################

# CG content is different from species to species
# If mix of bases changes along the read

def findGCByPos(reads):
    gc = [0] * 100
    totals = [0] * 100

    for read in reads:
        for i in range(len(read)):
            if read[i] == "C" or read[i] == "G":
                gc[i] += 1
            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc

gc = findGCByPos(seqs)
plt.plot(range(len(gc)), gc)

import collections
count = collections.Counter()
for seq in seqs:
    count.update(seq)
print(count)

############################################
# Sequencer give pieces to genomic puzzles #
############################################

'''
reads are far too short to tell story, it needs to be re-assembled into genome

unrelated humans have genomes that are 99.8-99.9% similarity

read alignment problem, where the read most closely matches to reference genome

problem arise when there is no reference genome -> assembly problem, de novo assembly
'''

####################################
# Read alignment and why it's hard #
####################################

'''
With help of reference genome, we can put the read back into a whole genome
1. with sequencing read, repeatedly going through for the best match
2. reads are few billions and reference genome can be x millions
'''

########################
# Naive exact matching #
########################

'''
Read alignment problem

1. Exact matching problem
At what offsets does pattern P occur within text T?
t = 'There would have been word'
t.find('word')
or
try every possible way to match 'word'

def naive(p,t): # naive algorithm
    occurences = []
    for i in range(len(t) - len(p) +1):
        match=True
        for j in range(len(p)):
            if t[t+j] != p[j]
            match=False
            break
        if match:
            occurences.append(i)
    return occurences

Let x = |P|, y = |T|
1. How many alignments are possible given x and y
-> y-x+1

2. What's the greatest # character comparisons possible?
-> x(y-x+1)

3. Worst case scenario on this naive algorithm is when there is a repeat matching every time.

4. What's the least # character comparision possible?
-> y-x+1
happens when first pattern of the text never occur
'''

#############################
# Matching artificial reads #
#############################

genome = readGenome('phix.fa')
def naive(p, t):
    '''
    p = read
    t = genome
    '''
    occurences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurences.append(i)
    return occurences
t = 'AGCTTAGATAGC'
p = 'AG'
print(naive(p,t))

# Generate artificial genome
import random
def generateReads(genome, numReads, readLen):
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome)-readLen) + 1
        reads.append(genome[start: start+readLen])
    return reads

reads = generateReads(genome, 100, 100)

numMatched = 0
for r in reads:
    matches = naive(r,genome)
    if len(matches) > 0:
        numMatched += 1
print("%d / %d reads matched exactly" % (numMatched, len(reads)))

#######################
# Matching read reads #
#######################

phix_reads, _ = readFASTQ('ERR266411_1.first1000.fastq')

numMatched = 0
n = 0
for r in phix_reads:
    r=r[:30]
    matches = naive(r,genome)
    matches.extend(naive(reverseComplement(r),genome))
    n+=1
    if len(matches) > 0:
        numMatched+=1
print("%d / %d reads matched genome" % (numMatched, n))

# it may be possible we are missing reverse complement
