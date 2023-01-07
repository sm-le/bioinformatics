# @sungml92
# Algorithm for assembly
# Programming homework 4

'''
1. In a practical, we saw the scs function (copied below along with overlap)
for finding the shortest common superstring of a set of strings.
'''
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

'''
It's possible for there to be multiple different shortest common superstrings for the same set of input strings.
Consider the input strings ABC, BCA, CAB. One shortest common superstring is ABCAB but another is BCABC and another is CABCA.

What is the length of the shortest common superstring of the following strings?
CCT, CTT, TGC, TGG, GAT, ATT
'''
strings = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]
len(scs(strings))

'''
2. How many different shortest common superstrings are there for the input strings given in the previous question?

Hint 1: You can modify the scs function to keep track of this.

Hint 2: You can look at these examples to double-check that your modified scs is working as expected.
'''

def multi_scs(ss):
    """ Returns multiple shortest common superstrings of given
        strings, which must be the same length """
    common_superstring = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        common_superstring.append(sup)
    shortest_length = min([len(i) for i in common_superstring])
    short_superstring = []
    for superstring in common_superstring:
        if len(superstring) == shortest_length:
            short_superstring.append(superstring)
    return short_superstring  # return shortest

multi_scs(strings)

'''
3. Download this FASTQ file containing synthetic sequencing reads from a mystery virus:

All the reads are the same length (100 bases) and are exact copies of substrings from the forward strand of the virus genome.
You don't have to worry about sequencing errors, ploidy, or reads coming from the reverse strand.

Assemble these reads using one of the approaches discussed, such as greedy shortest common superstring.
Since there are many reads, you might consider ways to make the algorithm faster,
such as the one discussed in the programming assignment in the previous module.

How many As are there in the full, assembled genome?

Hint: the virus genome you are assembling is exactly 15,894 bases long
'''

def pick_maximal_overlap(reads,k):
    readA, readB = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads,2):
        olen=overlap(a,b,min_length=k)
        if olen > best_olen:
            readA, readB = a, b
            best_olen = olen
    return readA, readB, best_olen

def greedy_scs(reads,k):
    readA, readB, olen = pick_maximal_overlap(reads,k)
    while olen>0:
        reads.remove(readA)
        reads.remove(readB)
        reads.append(readA+readB[olen:])
        readA, readB, olen = pick_maximal_overlap(reads,k)
    return ''.join(reads)

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

reads, _ = readFASTQ("ads1_week4_reads.fq")
synthetic_genome = greedy_scs(reads,8)

def countN(genome):
    count_dict = {}
    for nucleotide in genome:
        count_dict[nucleotide] = count_dict.get(nucleotide,0)+1
    return count_dict

'''
4. How many Ts are there in the full, assembled genome from the previous question?
'''
# -*- coding: utf-8 -*-
"""
imported from zhihuixie's work
"""
import itertools

class deNovoAssambly ():
    """
    This class is used to assembly genomic sequence from sequencing reads.
    """
    def __init__(self, filename = False, text = False):
        # initiate paramater
        self.filename = filename
        self.text = text
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

    def phraseReads(self, reads, k_mer):
        """
        construct the prefix and suffix of a read to a dictornary with read as
        key and pre-,suffix as values
        """
        reads_dict = {}
        for read in reads:
            for i in range(len(read) - k_mer + 1):
                substring = read[i:i+k_mer]
                if substring not in reads_dict:
                    reads_dict[substring] = set([read])
                else:
                    reads_dict[substring].add(read)
        return reads_dict

    def bestOverlapFast(self, reads, k_mer):
        """
        return the reads pair and maxium length of the overlapped pair in reads
        """
        reads_dict = self.phraseReads(reads, k_mer)
        read_a, read_b = None, None
        best_overlap = 0
        for read1 in reads:
            k_mer_string = read1[len(read1) - k_mer:]
            if k_mer_string in reads_dict:
                reads_set = reads_dict[k_mer_string]
                for read2 in reads_set:
                    if read1 != read2: #skip self comparison
                        offset = self.overlap(read1, read2, k_mer)
                        if offset > best_overlap: # skip non-overlapped pairs
                           best_overlap = offset
                           read_a, read_b = read1, read2
        return read_a, read_b, best_overlap

    def shortComSupstr(self):
        """ Returns shortest common superstring and its length of given
        strings, which must be the same length """
        text = self.text
        sup_strs = []
        shortest_length = None
        for string_pair in itertools.permutations(text):
            sup_str = string_pair[0]
            for i in range(len(text) - 1):
                overlap_length = self.overlap(string_pair[i], string_pair[i+1], 1)
                sup_str += string_pair[i+1][overlap_length:]
            sup_strs.append(sup_str)
            if shortest_length is None or shortest_length >= len(sup_str):
                shortest_length = len(sup_str)
        shortest_sup = [s for s in sup_strs if len(s) == shortest_length]
        return shortest_length, shortest_sup

    def de_bruijnGraph(self, k):
        """
        build a de bruijn graph for a give sequence
        """
        sequences, _ = self.readFastq()
        print("read file completed...")
        graph = {}
        for sequence in sequences:
            for i in range(len(sequence) - k + 1):
                node = sequence[i:i+k-1]
                if node not in graph:
                    graph[node] = [sequence[i+1:i+k]]
                else:
                    graph[node].append(sequence[i+1:i+k])
        return graph

    def greedyScs(self, reads, k_mer):
        """
        use greedy algorithm to find the maxi overlap sequence
        """
        read1, read2, best_overlap = self.bestOverlapFast(reads, k_mer)
        while best_overlap > 0:
            reads.remove(read1)
            reads.remove(read2)
            reads.append(read1 + read2[best_overlap:])
            read1, read2, best_overlap = self.bestOverlapFast(reads, k_mer)

        return "".join(reads)

    def output(self, sequence):
        """
        output a given sequence to a txt file
        """
        with open ("assembly sequence of " + self.filename + ".txt", "w") as f:
            f.write(sequence)

    def countNt(self,  sequence):
        """count number of each base in a given sequence
        """
        dict_count = {"A":0, "C":0, "G":0, "T":0, "N":0}
        for nt in sequence:
            if nt in dict_count:
                dict_count[nt] = dict_count.get(nt) + 1
            else:
                dict_count[nt] = 1
        return dict_count


if __name__ == "__main__":
    # test case1
    text = ['ABC', 'BCA', 'CAB']
    # test case2
    #text = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
    # questions 1 and 2
    #text = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]
    #assembly = deNovoAssambly (text = text)
    #print assembly.shortComSupstr()
    filename = "ads1_week4_reads.fq"
    assembly = deNovoAssambly (filename = filename)
    reads, _ = assembly.readFastq()
    sequence = assembly.greedyScs(reads, 8)
    counts = assembly.countNt(sequence)
    assembly.output(sequence)
    print(counts)
    """
    from sequencer import *
    dg = DeBruijn(filename)
    path = dg.sequence()
    print path, len(path), 15894
    """
