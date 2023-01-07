'''
In a practical, we saw Python code implementing the Boyer-Moore algorithm.
Some of the code is for preprocessing the pattern P into the tables needed to execute the bad character and good suffix rules —
we did not discuss that code. But we did discuss the code that performs the algorithm given those tables:
'''
import bm_preproc

def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences
'''
Measuring Boyer-Moore's benefit. First, download the Python module for Boyer-Moore preprocessing:

This module provides the BoyerMoore class, which encapsulates the preprocessing info used by the boyer_moore function above.
Second, download the provided excerpt of human chromosome 1:


Third, implement versions of the naive exact matching and Boyer-Moore algorithms that additionally count and return
(a) the number of character comparisons performed and
(b) the number of alignments tried.
Roughly speaking, these measure how much work the two different algorithms are doing.
'''

def naive_with_counts(p, t):
    occurrences = []
    num_aln = 0
    num_char_comp = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        num_aln += 1
        match = True
        for j in range(len(p)):  # loop over characters
            num_char_comp += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_aln, num_char_comp

'''
For a few examples to help you test if your enhanced versions of
the naive exact matching and Boyer-Moore algorithms are working properly, see these notebooks:

# Implement naive_with_counts by extending naive function
from naive_with_counts import naive_with_counts

Example 1
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)
([40], 41, 46)

Example 2
p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)
([0, 19], 20, 35)
'''
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

'''
# Implement boyer_moore_with_counts by extending boyer_moore function
from bm_with_counts import boyer_moore_with_counts
from bm_preproc import BoyerMoore
'''

def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    num_aln = 0
    num_char_comp = 0
    while i < len(t) - len(p) + 1:
        num_aln += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            num_char_comp += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_aln, num_char_comp

'''
Example1
p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)
([40], 12, 15)

Example 2
p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)
([0, 19], 5, 18)
'''
from bm_preproc import BoyerMoore
p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)


'''
1. How many alignments does the naive exact matching algorithm try when matching the string
GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)

2. How many character comparisons does the naive exact matching algorithm try when matching the string
GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)
'''

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = readGenome('chr1.GRCh38.excerpt.fasta')

naive_with_counts(p,t)

'''
3. How many alignments does Boyer-Moore try when matching the string
GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)
'''
p_bm = BoyerMoore(p, lowercase_alphabet.upper())
boyer_moore_with_counts(p,p_bm,t)

'''
4. Index-assisted approximate matching. In practicals, we built a Python class called Index
implementing an ordered-list version of the k-mer index. The Index class is copied below.

We also implemented the pigeonhole principle using Boyer-Moore as our exact matching algorithm.

Implement the pigeonhole principle using Index to find exact matches for the partitions.
Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions). We will use an 8-mer index.

Download the Python module for building a k-mer index.
https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py

Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches.
Insertions and deletions are not allowed. Don't consider any reverse complements.

How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence,
occur with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse complements here.)

Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count a match more than once.

Hint 2: You can check your work by comparing the output of your new function to that of the naive_2mm function implemented in the previous module.
'''

import kmer_index
import bisect

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
t = readGenome('chr1.GRCh38.excerpt.fasta')
index = kmer_index.Index(t,8)

def pigeonhole_index_match(p,t,n,index): # pigeonhole principle divide p by n+1
    segment_length = round(len(p)/(n+1))
    all_matches = set()
    index_hits = []
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = index.query(p[start:end])
        for m in matches:
            index_hits.append(m)
            if m < start or m-start+len(p) > len(t): # location stay in start to end
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches+=1
                    if mismatches>n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches+=1
                    if mismatches>n:
                        break
            if mismatches <= n:
                all_matches.add(m-start)
    return list(all_matches), index_hits

list_of_matches, list_of_hits = pigeonhole_index_match(p,t,2,index)

print(sorted(list_of_hits))
print(len(list_of_matches))
'''
5. Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of
GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse complements.)

Hint: You should be able to use the boyer_moore function (or the slower naive function) to double-check your answer.
'''

print(len(list_of_hits))


'''
Let's examine whether there is a benefit to using an index built using subsequences of T rather than substrings,
as we discussed in the "Variations on k-mer indexes" video. We'll consider subsequences involving every N characters.
For example, if we split ATATAT into two substring partitions, we would get partitions ATA (the first half)
and TAT (second half). But if we split ATATAT into two subsequences by taking every other character,
we would get AAA (first, third and fifth characters) and TTT (second, fourth and sixth).

Another way to visualize this is using numbers to show how each character of P is allocated to a partition.
Splitting a length-6 pattern into two substrings could be represented as 111222,
and splitting into two subsequences of every other character could be represented as 121212

The following class SubseqIndex is a more general implementation of Index that additionally handles subsequences.
It only considers subsequences that take every Nth character:
'''

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

subIndex = SubseqIndex(t,8,3)

def pigeonhole_subindex_match(p,t,n,subIndex): # pigeonhole principle divide p by n+1
    segment_length = round(len(p)/(n+1))
    all_matches = set()
    index_hits = []
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = subIndex.query(p[start:end])
        index_hits.append(matches)
        for m in matches:
            if m < start or m-start+len(p) > len(t): # location stay in start to end
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches+=1
                    if mismatches>n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches+=1
                    if mismatches>n:
                        break
            if mismatches <= n:
                all_matches.add(m-start)
    return list(all_matches), index_hits

list_of_submatches, no_subHits = pigeonhole_subindex_match(p,t,2,subIndex)
print(len(list_of_submatches))

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = pigeonhole_subindex_match(p, t, 2,subseq_ind)
print(occurrences)
print(num_index_hits)






'''


'''


import bisect

class findPatternV2 ():
    """
    This class finds the occurance and position of a given pattern in a given
    genomic sequence in a file.
    """
    def __init__(self, pattern, filename = False):
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

    def naiveMatch (self, numberOfMismatch, text = False):
        """
        this is naive match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        if text != False: #for test cases
            genome = text
        else:
            genome = self.readGenome()
        pattern = self.pattern
        occurances = []
        alignments = 0
        comparisons = 0
        for i in range(len(genome) - len(pattern) + 1):
            match = True
            counter = 0
            for j in range(len(pattern)):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    counter += 1
                if counter > numberOfMismatch:
                    match = False
                    break
            if match:
                occurances.append(i)
            alignments += 1
        return occurances, alignments, comparisons

    def boyerMoore (self, numberOfMismatch, bm, text = False):
        """
        this is naive match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        i = 0
        if text != False: #for test cases
            genome = text
        else:
            genome = self.readGenome()
        pattern = self.pattern
        occurances = []
        alignments = 0
        comparisons = 0
        while i < len(genome) - len(pattern) + 1:
            shift = 1
            match = True
            for j in range(len(pattern) - 1, -1, -1):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    badCharacterSkip = bm.bad_character_rule(j, genome[i+j])
                    goodSuffixSkip = bm.good_suffix_rule(j)
                    shift = max(shift, badCharacterSkip, goodSuffixSkip)
                    match = False
                    break
            if match:
                occurances.append(i)
                goodSuffixSkip = bm.match_skip()
                shift = max(shift, goodSuffixSkip)
            i += shift
            alignments += 1
        return occurances, alignments, comparisons

    def matchedIndex(self, index, k_mer, pattern, isSubseqIndex):
        """
        find number of hits, occurances and time of occurance for a given pattern
        using string index in genome
        """
        genome = self.readGenome()
        occurances_match = []
        hit_index = []
        occurance_genome = []
        counter = 0
        if not isSubseqIndex:
            length = len(pattern)-k_mer + 1
        else:
            length = isSubseqIndex
        for i in range(length): # loop over to generate kmers
            if not isSubseqIndex:
                pattern_q = pattern[i:i+k_mer]
            else:
                pattern_q = pattern[i:]
            hits = index.query(pattern_q) # query each kmer
            for hit in hits:
                counter += 1 #count total number of hits
                text = genome[hit-i : hit+len(pattern) -i]
                if hit-i not in hit_index: #avoid duplicated counts
                    hit_index.append(hit-i)
                    occurance, _, _ = self.naiveMatch(2,text)
                    occurances_match.extend(occurance)
                if len(occurance) != 0 and hit-i not in occurance_genome:
                    occurance_genome.append(hit-i)
        return occurance_genome, len(occurances_match), counter

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    def genome_index(self):
         return self.index

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

if __name__ == "__main__":
    from bm_preproc import BoyerMoore
    #Questions 1-3
    filename='chr1.GRCh38.excerpt.fasta'
    #filename = ("../data/chr1.GRCh38.excerpt.fasta")
    #Q1: How many alignments does the naive exact matching algorithm try when
    #matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
    #(derived from human Alu sequences) to the excerpt of human chromosome 1?
    #(Don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    patterns = findPatternV2(pattern, filename)
    print("Q1: The alignments for naive match algorithm is", patterns.naiveMatch(0)[1])
    patterns = findPatternV2("GGCGCGGTGGCTCACGCCTGTAAT", filename)

    #Q2: How many character comparisons does the naive exact matching algorithm
    #try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
    #(derived from human Alu sequences) to the excerpt of human chromosome 1?
    #(Don't consider reverse complements.)
    print("Q2: The characters comparisons for naive match algorithm is %d" % patterns.naiveMatch(0)[2])

    #How many alignments does Boyer-Moore try when matching the string
    #GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
    #(derived from human Alu sequences) to the excerpt of human chromosome 1?
    #(Don't consider reverse complements.)
    print("Q3: The alignments for Boyer-Moore algorithm is %d" % patterns.boyerMoore(0, \
          BoyerMoore(pattern, "ACGT"))[1])

    #Q4: How many times does the string GGCGCGGTGGCTCACGCCTGTAAT,
    #which is derived from a human Alu sequence, occur with up to 2
    #substitutions in the excerpt of human chromosome 1?
    #(Don't consider reverse complements here.)
    k_mer = 8
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    genome = patterns.readGenome()
    index = Index(genome, k_mer)
    occurances, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, \
                                               k_mer, pattern, isSubseqIndex = False)
    print("Q4: Within 2 mismatchs, the string occurs %d times\n" %numberOfOccurs)

    #Q5:Using the instructions given in Question 4, how many total index hits
    #are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with
    #up to 2 substitutions in the excerpt of human chromosome 1?
    print("Q5: Within 2 mismatchs, the total index hits are %d \n" %numberOfhits)

    #Q6: When using this function, how many total index hits are there when
    #searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the
    #excerpt of human chromosome 1? (Again, don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    k_mer = 8
    vial = 3
    index = SubseqIndex(t, k_mer, vial)
    occurances, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, \
                                               k_mer, pattern, isSubseqIndex = vial)

    print("Q6: Within 2 mismatchs, the hits are %d\n" %numberOfhits)
