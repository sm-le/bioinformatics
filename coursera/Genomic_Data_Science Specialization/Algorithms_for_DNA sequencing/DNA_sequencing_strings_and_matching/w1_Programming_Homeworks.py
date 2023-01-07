# In lecture and in a practical, we saw an implementation of the naive exact matching algorithm:

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# ...and we saw a function that takes a DNA string and returns its reverse complement:

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

# ...and we saw a function that parses a DNA reference genome from a file in the FASTA format.
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
# ...and we saw a function that parses the read and quality strings from a FASTQ file containing sequencing reads.
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
# First, implement a version of the naive exact matching algorithm that is strand-aware.
# That is, instead of looking only for occurrences of P in T,
# additionally look for occurrences of thereverse complement of P in T. If P is ACT,
# your function should find occurrences of both ACTand its reverse complement AGT in T.

# If P and its reverse complement are identical (e.g. AACGTT),
# then a given match offset should be reported only once.
# So if your new function is called naive_with_rc,
# then the old naivefunction and your new naive_with_rc function should return the same results when P equals its reverse complement.

def naive_with_rc(p, t):
    rev_p = reverseComplement(p)
    occurrences = []

    for i in range(len(t)-len(p)+1):
        match = True

        for j in range(len(p)):
            if (t[i+j] != p[j]) and (t[i+j] != rev_p[j]):
                match = False
                break

        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naive_with_rc(p, t):
    rev_p = reverseComplement(p)
    occurrences = []

    for i in range(len(t)-len(p)+1):
        match = True

        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)

        if p != rev_p:
            match = True
            for j in range(len(p)):
                if t[i+j] != rev_p[j]:
                    match = False
                    break

        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences




################
# testing code #
################

'''
Example 1:

Expected result: [10,23]
'''
p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

'''
Example 2:

[10,24]
'''
p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

'''
Example 3

offset of leftmost occurrence: 62
# occurrences: 60
'''

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

phix_genome = readGenome('phix.fa')
occurrences = naive_with_rc('ATTA', phix_genome)
print(occurrences)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))

'''
1. How many times does AGGT or its reverse complement ACCT occur in the lambda virus genome?
E.g. if AGGT occurs 10 times and ACCT occurs 12 times, you should report 22.
p
'''

lambda_genome = readGenome('lambda_virus.fa')
occurences = naive_with_rc('AGGT',lambda_genome)
print(len(occurences))

'''
2. How many times does TTAA or its reverse complement occur in the lambda virus genome?

Hint: TTAA and its reverse complement are equal, so remember not to double count.
p
'''

occurences = naive_with_rc('TTAA',lambda_genome)
print(occurences)
print(len(occurences))

'''
3. What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?
E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement ACTTAGT is at offset 29, then report 29.
p
'''
occurences = naive_with_rc('ACTAAGT',lambda_genome)
print(min(occurences))


'''
4. What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?
p

'''
occurences = naive_with_rc('AGTCGA',lambda_genome)
print(min(occurences))


'''
5. As we will discuss, sometimes we would like to find approximate matches for P in T. That is, we want to find occurrences with one or more differences.

For Questions 5 and 6, make a new version of the naive function called naive_2mm that allows up to 2 mismatches per occurrence. Unlike for the previous questions, do not consider the reverse complement here. We're looking for approximate matches for P itself, not its reverse complement.

￼

For example, ACTTTA occurs twice in ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, and once at offset 4 with 1 mismatch. So naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT') should return the list [0,4].

Hint: See this notebook for a few examples you can use to test your naive_2mm function.

How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?
p
'''

def naive_2mm(p,t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        count = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                count +=1
        if count > 2:
            match = False
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


occurrences = naive_2mm('TTCAAGCC',lambda_genome)
print(len(occurrences))

'''
6. What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?
p
'''
occurrences = naive_2mm('AGGAGGTT',lambda_genome)
print(occurrences)
print(min(occurrences))

'''
7. Finally, download and parse the provided FASTQ file containing real DNA sequencing reads derived from a human:

https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq

Note that the file has many reads in it and you should examine all of them together when answering this question. The reads are taken from this study:

Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011). Accurate

and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505.

This dataset has something wrong with it; one of the sequencing cycles is poor quality.

Report which sequencing cycle has the problem. Remember that a sequencing cycle corresponds to a particular offset in all the reads. For example, if the leftmost read position seems to have a problem consistently across reads, report 0. If the fourth position from the left has the problem, report 3. Do whatever analysis you think is needed to identify the bad cycle. It might help to review the "Analyzing reads by position" video.
f
'''

seqs, quals = readFastq("ERR037900_1.first1000.fastq")
print(quals)
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

def findPosBQual(quals):
    occurrences = []
    for qual in quals:
        for i in range(len(qual)):
            qc = phred33ToQ(qual[i])
            if qc < 3:
                occurrences.append(i)

    return occurrences

f = findPosBQual(quals)
plt.hist(f,bins=100)

def most_frequent(List):
    return max(set(List), key = List.count)
most_frequent(f)




# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:42:50 2015
This project includes two classes: findPattern and checkQuality. The findPattern
class finds the occurance and position of a given pattern in a given genomic
sequence. The checkQuality class exams quality of sequencing for each cycle
@author: zhihuixie
"""


class findPattern ():
    """
    This class finds the occurance and position of a given pattern in a given
    genomic sequence in a file.
    """
    def __init__(self, pattern, filename):
        # initiate parameters
        self.pattern = pattern
        self.filename = filename

    def readGenome (self):
        """
        read genomic DNA sequence to a string
        """
        genome = ""
        with open (self.filename, "r") as f:
            for line in f:
                # skip header
                if not line[0] == ">":
                    genome += line.rstrip()
            f.close()
        return genome


    def reverseComplement (self):
        """
        generate reverse complement sequence for a given dna sequence
        """
        complement = {"A": "T", "C": "G", "T": "A", "G": "C"}
        revComPattern = "" # reversed compliment pattern
        for nt in self.pattern:
            revComPattern = complement[nt] + revComPattern

        return revComPattern

    def match(self, string1, string2, numOfMismatch):
        """
        return True or False for matching results of two strings under the offset
        of numOfMismatch
        """
        counter = 0
        if len(string1) != len(string2):
            return False
        # loop over string to compare character
        for i in range(len(string1)):
            if string1[i] != string2[i]:
                counter += 1
        if counter > numOfMismatch:
            return False
        return True

    def patternIdentifier (self, numOfMismatch):
        """
        find positions of a given pattern and the reversed complement pattern
        in a given genome
        """
        patternLength = len(self.pattern)
        genome = self.readGenome()
        revComPattern = self.reverseComplement()
        occurances = []

        for i in range (patternLength): # loop over pattern index
            # loop over genome patterns
            for j in range (i, len(genome), patternLength):
                genomeMotif = genome[j: j+patternLength]
                # compare genomic motif and patterns
                if numOfMismatch == 0:
                    if (self.match(genomeMotif, self.pattern, 0) or \
                        self.match(genomeMotif, revComPattern, 0))\
                       and j not in occurances: # avoid duplicate records
                        occurances.append(j)
                else:
                    if self.match(genomeMotif, self.pattern, numOfMismatch)\
                       and j not in occurances: # avoid duplicate records:
                        occurances.append(j)
        return occurances

class checkQuality ():
    """
    The checkQuality class exams quality of sequencing for each cycle
    """
    def __init__ (self, filename):
        self.filename = filename

    def readFastq(self):
        """
        read dna sequence and quality base from a fastq sequencing file to lists
        """
        with open (self.filename, "r") as f:
             sequences = []
             qualities = []
             while True:
                 f.readline() # skip name line
                 seq = f.readline().rstrip() # read sequence line
                 f.readline() # skip strand line
                 qual = f.readline().rstrip() # read quality line
                 if len(seq) == 0: #finish read
                     break
                 # add seqence and quality information to list
                 sequences.append(seq)
                 qualities.append(qual)
             f.close()
        return sequences, qualities

    def phre33ToQ (self, qualString):
        """
        transform quality string to quality base score
        """
        qScore = []
        for qual in qualString:
            qScore.append(ord(qual) - 33)
        return qScore
    def findPoorQuality(self):
        """
        find the index of poorest q score in each sequencing
        """
        _, qualities = self.readFastq()
        lowestQScoreIndex = []
        for qualString in qualities:
            qScore = self.phre33ToQ(qualString)
            lowestQScoreIndex.append(qScore.index(min(qScore)))
        return lowestQScoreIndex

    def countPoorQuality(self):
        """
        count number of poorest q score in each cycle
        """
        import collections
        return collections.Counter(self.findPoorQuality())

    def plotHist(self):
        """
         show the distribution of poorest q score
        """
        import matplotlib.pyplot as plt
        data = self.countPoorQuality()
        plt.bar(data.keys(), data.values())
        plt.show()


if __name__ == "__main__":
    #Test
    filename = "../data/phix.fa"
    pattern = "ATTA"
    patterns = findPattern(pattern,filename)
    print "Test dataset results - Occurances and leftmost offset: "
    print len(patterns.patternIdentifier(0)), min(patterns.patternIdentifier(0)), "\n"

    #Questions 1-6
    filename1 = "../data/lambda_virus.fa"
    #Q1: How many times does AGGT or its reverse complement (ACCT) occur in the
    #lambda virus genome? E.g. if AGGT occurs 10 times and ACCT occurs 12 times,
    #you should report 22.
    pattern = "AGGT"
    patterns = findPattern(pattern,filename1)
    print "Q1: The 'AGGT' or 'ACCT' occurs %d times \n" \
           %len(patterns.patternIdentifier(0))

    #Q2: How many times does TTAA or its reverse complement occur in the lambda
    #virus genome? Hint: TTAA and its reverse complement are equal, so remember
    #not to double count.
    pattern = "TTAA"
    patterns = findPattern(pattern,filename1)
    print "Q2: The 'TTAA' occurs %d times \n" \
           %len(patterns.patternIdentifier(0))

    #Q3: What is the offset of the leftmost occurrence of ACTAAGT or its reverse
    #complement in the Lambda virus genome? E.g. if the leftmost occurrence of
    #ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse
    #complement ACTTAGT is at offset 29, then report 29.
    pattern = "ACTAAGT"
    patterns = findPattern(pattern,filename1)
    print "Q3: The offset of the leftmost occurrence of ACTAAGT is %d \n" \
           %min(patterns.patternIdentifier(0))

    #Q4: What is the offset of the leftmost occurrence of AGTCGA or its reverse
    #complement in the Lambda virus genome?
    pattern = "AGTCGA"
    patterns = findPattern(pattern,filename1)
    print "Q4: The offset of the leftmost occurrence of AGTCGA is %d \n" \
           %min(patterns.patternIdentifier(0))

    #Q5: How many times does TTCAAGCC occur in the Lambda virus genome when
    #allowing up to 2 mismatches?
    pattern = "TTCAAGCC"
    patterns = findPattern(pattern,filename1)
    print "Q5: The 'TTCAAGCC' occurs %d times with mismatch less than 2 \n" \
           %len(patterns.patternIdentifier(2))

    #Q6: What is the offset of the leftmost occurrence of AGGAGGTT in the
    #Lambda virus genome when allowing up to 2 mismatches?
    pattern = "AGGAGGTT"
    patterns = findPattern(pattern,filename1)
    print "Q6: The offset of the leftmost occurrence of AGGAGGTT with mismatch less than 2 is %d \n" \
          %min(patterns.patternIdentifier(2))

    #Q7: Report which sequencing cycle has the problem.
    filename2 = "../data/ERR037900_1.first1000.fastq"
    qualities = checkQuality(filename2)
    counters = qualities.countPoorQuality()
    print "Q7: The cycle has most frequent poor quality is %d" \
          %max(counters, key = lambda x: counters[x])
    
