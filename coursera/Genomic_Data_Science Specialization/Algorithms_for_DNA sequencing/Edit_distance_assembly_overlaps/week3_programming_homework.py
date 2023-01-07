'''
We saw how to adapt dynamic programming to find approximate occurrences of a pattern in a text. Recall that:

Rows of the dynamic programming matrix are labeled with bases from P and columns with bases from T
Elements in the first row are set to 0
Elements in the first column are set to 0, 1, 2, ..., as for edit distance
Other elements are set in the same way as elements of a standard edit distance matrix
The minimal value in the bottom row is the edit distance of the closest match between P and T

First, download the provided excerpt of human chromosome 1

Second, parse it using the readGenome function we wrote before.

Third, adapt the editDistance function we saw in practical (copied below) to answer questions 1 and 2 below.

Your function should take arguments p (pattern), t (text) and should return
the edit distance of the match between P and T with the fewest edits.


Hint: In the "A new solution to approximate matching" video we saw that the best approximate match of P = GCGTATGC within
T = TATTGGCTATACGGTT had 2 edits. You can use this and other small examples to double-check that your function is working.
'''
def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1][:])

editDistance('GCGTATGC','TATTGGCTATACGGTT')

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

'''
1. What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1? (Don't consider reverse complements.)
'''
P = 'GCTGATCGATCGTACG'
T = readGenome('chr1.GRCh38.excerpt.fasta')
editDistance(P,T)


'''
2. What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1? (Don't consider reverse complements.)
'''

P = 'GATTTACCAGATTGAG'
T = readGenome('chr1.GRCh38.excerpt.fasta')
editDistance(P,T)


'''
3. In a practical, we saw a function for finding the longest exact overlap (suffix/prefix match) between two strings. The function is copied below.
'''

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

'''
Say we are concerned only with overlaps that (a) are exact matches (no differences allowed), and (b) are at least k bases long.
To make an overlap graph, we could call overlap(a,b,min_length=k) on every possible pair of reads from the dataset.
Unfortunately, that will be very slow!

Consider this: Say we are using k=6, and we have a read a whose length-6 suffix is GTCCTA.
Say GTCCTA does not occur in any other read in the dataset.
In other words, the 6-mer GTCCTA occurs at the end of read a and nowhere else.
It follows that a's suffix cannot possibly overlap the prefix of any other read by 6 or more characters.

Put another way, if we want to find the overlaps involving a suffix of read a and a prefix of some other read,
we can ignore any reads that don't contain the length-k suffix of a. This is good news because it can save us a lot of work!

Here is a suggestion for how to implement this idea. You don't have to do it this way, but this might help you.
Let every k-mer in the dataset have an associated Python set object, which starts out empty.
We use a Python dictionary to associate each k-mer with its corresponding set.
(1) For every k-mer in a read, we add the read to the set object corresponding to that k-mer.
If our read is GATTA and k=3, we would add GATTA to the set objects for GAT,
ATT and TTA. We do this for every read so that, at the end, each set contains all reads containing the corresponding k-mer.
(2) Now, for each read a, we find all overlaps involving a suffix of a. To do this, we take a's length-k suffix,
find all reads containing that k-mer (obtained from the corresponding set) and call overlap(a,b,min_length=k) for each.

The most important point is that we do not call overlap(a,b,min_length=k) if b does not contain the length-k suffix of a.

Download and parse the read sequences from the provided Phi-X FASTQ file. We'll just use their base sequences,
so you can ignore read names and base qualities. Also, no two reads in the FASTQ have the same sequence of bases. This makes things simpler.

https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.for_asm.fastq

Next, find all pairs of reads with an exact suffix/prefix match of length at least 30. Don't overlap a read with itself;
if a read has a suffix/prefix match to itself, ignore that match. Ignore reverse complements.

Hint 1: Your function should not take much more than 15 seconds to run on this 10,000-read dataset, and maybe much less than that. (Our solution takes about 3 seconds.)
If your function is much slower, there is a problem somewhere.
Hint 2: Remember not to overlap a read with itself. If you do, your answers will be too high.
Hint 3: You can test your implementation by making up small examples, then checking that (a) your implementation runs quickly,
and (b) you get the same answer as if you had simply called overlap(a,b,min_length=k) on every pair of reads.
We also have provided a couple examples you can check against.

Picture the overlap graph corresponding to the overlaps just calculated. How many edges are in the graph? In other words, how many distinct pairs of reads overlap?
'''
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

seqs, _ = readFASTQ('ERR266411_1.for_asm.fastq')
kmer_dict = {}
kmer=30
for seq in seqs:
    for pos in range(len(seq)-kmer+1):
        substring = seq[pos:pos+kmer]
        if substring not in kmer_dict:
            kmer_dict[substring] = set([seq])
        else:
            kmer_dict[substring].add(seq)
kmer_dict

overlap_graph = {}
for seq in seqs:
    suffix_kmer = seq[len(seq)-kmer:]
    if suffix_kmer in kmer_dict:
        edges = set([])
        reads_set = kmer_dict[suffix_kmer]
        for read in reads_set:
            if seq != read:
                overlapped = overlap(seq,read,kmer)
                if overlapped > 0:
                    edges.add(read)
                    overlap_graph[seq] = edges

len(overlap_graph)
sum([len(i) for i in overlap_graph.values()])

def overlapGraph(self, k_mer):
    """
    construct graph with key as a read (node) and values as all other
    reads overlapped with the previous read (node)
    """
    reads_dict = self.phraseReads(k_mer)
    reads, _= self.readFastq()
    graph = {}
    for read1 in reads:
        k_mer_string = read1[len(read1) - k_mer:]
        if k_mer_string in reads_dict:
            edges = set([])
            reads_set = reads_dict[k_mer_string]
            for read2 in reads_set:
                if read1 != read2: #skip self comparison
                    offset = self.overlap(read1, read2, k_mer)
                    if offset > 0: # skip non-overlapped pairs
                        edges.add(read2) #add overlapped reads to be values
                        graph[read1] = edges
    return graph

'''
4. Picture the overlap graph corresponding to the overlaps computed for the previous question.
How many nodes in this graph have at least one outgoing edge? (In other words, how many reads have a suffix involved in an overlap?)
'''







# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 19:55:08 2015
This project includes findPatternV3 class. The findPatternV3
class finds the edit distance of a given pattern in a given genomic
sequence and constructs overlap graphs.
@author: zhihuixie
"""

from itertools import permutations

class findPatternV3 ():
    """
    This class finds the edit distance of a given pattern in a given genomic
    sequence and constructs overlap graphs.
    """
    def __init__(self, filename, pattern = False):
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

    def editDistance(self):
        """
        Implement dynamic algorithm to calculate edit distance between a given
        pattern and a given genome
        """
        pattern = self.pattern
        genome = self.readGenome()
        pattern_length = len(pattern) + 1
        genome_length = len(genome) + 1
        #generate matrix
        matrix = [[0]*genome_length for i in range(pattern_length)]
        # initiate the first column
        for i in range(pattern_length):
            matrix[i][0] = i
        for i in range(1, pattern_length):
            for j in range(1, genome_length):
                dist_hor = matrix[i][j-1] + 1
                dist_vel = matrix[i-1][j] + 1
                dist_diag = matrix[i-1][j-1] + 1 if pattern[i-1] != genome[j-1]\
                            else matrix[i-1][j-1]
                matrix[i][j] = min(dist_hor, dist_vel, dist_diag)
        return min(matrix[-1])

    def phraseReads(self, k_mer):
        """
        construct the prefix and suffix of a read to a dictornary with read as
        key and pre-,suffix as values
        """
        reads, _ = self.readFastq()
        reads_dict = {}
        for read in reads:
            for i in range(len(read) - k_mer + 1):
                substring = read[i:i+k_mer]
                if substring not in reads_dict:
                    reads_dict[substring] = set([read])
                else:
                    reads_dict[substring].add(read)
        return reads_dict
    def overlap(self, read1, read2, k_mer):
        """
        find overlaped leftmost offset
        """
        start = 0
        while True:
            start = read1.find(read2[:k_mer], start)
            if start == -1:
                return 0 # without overlap
            if read2.startswith(read1[start:]):
                return len(read1) - start
            start += 1

    def overlapGraph(self, k_mer):
        """
        construct graph with key as a read (node) and values as all other
        reads overlapped with the previous read (node)
        """
        reads_dict = self.phraseReads(k_mer)
        reads, _= self.readFastq()
        graph = {}
        for read1 in reads:
            k_mer_string = read1[len(read1) - k_mer:]
            if k_mer_string in reads_dict:
                edges = set([])
                reads_set = reads_dict[k_mer_string]
                for read2 in reads_set:
                    if read1 != read2: #skip self comparison
                        offset = self.overlap(read1, read2, k_mer)
                        if offset > 0: # skip non-overlapped pairs
                            edges.add(read2) #add overlapped reads to be values
                            graph[read1] = edges
        return graph

    def naive_overlap_map(self, k_mer):
        """
        construct graph with key as a pair of reads with overlap and values as
        the leftmost offset of the overlap
        """
        graph = {}
        reads, _ = self.readFastq()
        for read1, read2 in permutations(reads, 2):
            #skip non-overlapped reads
            if read1[len(read1) - k_mer:] in read2:
                offset = self.overlap(read1, read2, k_mer)
                    # check if reads[i] overlapped with reads[j]
                if offset != 0:
                    graph[(read1,read2)] = offset
        return graph

if __name__ == "__main__":

    #Q1: What is the edit distance of the best match between pattern
    #GCTGATCGATCGTACG and the excerpt of human chromosome 1?
    #(Don't consider reverse complements.)
    pattern = "GCTGATCGATCGTACG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3 (filename, pattern)
    edit_dist = patterns.editDistance()
    print "Q1: the edit distance of the best match between pattern and the genome is %d\n"\
           %edit_dist

    #Q2: What is the edit distance of the best match between pattern
    #GATTTACCAGATTGAG and the excerpt of human chromosome 1?
    #(Don't consider reverse complements.)
    pattern = "GATTTACCAGATTGAG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3 (filename, pattern)
    edit_dist = patterns.editDistance()
    print "Q2: the edit distance of the best match between pattern and the genome is %d\n"\
           %edit_dist

    #Q3: Picture the overlap graph corresponding to the overlaps just calculated.
    # How many edges are in the graph? In other words, how many distinct pairs
    # of reads overlap?
    #Q4: Picture the overlap graph corresponding to the overlaps computed for
    #the previous question. How many nodes in this graph have at least one
    #outgoing edge? (In other words, how many reads have a suffix involved in
    #an overlap?)
    import time
    t1 = time.time()
    filename = "../data/ERR266411_1.for_asm.fastq"
    patterns = findPatternV3 (filename)
    k_mer = 30
    graph = patterns.naive_overlap_map(k_mer)
    t2 = time.time()
    print "Running time for naive overlap mapping: %d sec\n"%(t2 - t1)

    reads = patterns.phraseReads(k_mer)
    t3 = time.time()
    print "Running time for phrase reads: %d sec\n"%(t3 - t2)
    graph = patterns.overlapGraph(k_mer)
    t4 = time.time()
    print "Running time for optimized algorithm: %d sec\n"%(t4 - t2)
    numberOfNodes = len(graph)
    numberOfEdges = sum([len(edges) for edges in graph.values()])
    print "Q3: the total edges are %d\n"%numberOfEdges
    print "Q4: the total nodes are %d"%numberOfNodes
