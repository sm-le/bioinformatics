# @sungml92
# Algorithm for DNA sequencing
# Week 2: Preprocessing, indexing, and approximate matching

################
# Introduction #
################

'''
Objective
1. work with a range of ideas that underlying fast pattern matching algorithm and both exact and approximate matching algorithm
2. Boyer-Moore
3. pigeon hole principle
'''

######################
# Boyer-Moore basics #
######################

'''
Boyer-Moore is similar to naive exact matching
- It will try characeter comparison in an alignment
- But it will skip alignment that does not need to examine


Exact matching: better naive algorithm

P: word
T: There would have been a time for such a word
         word

u doesn't occur in P, so we can skip next two alignments

For example,
T: There would have been a time for such a word
         word
          word - skip
           word - skip

Boyer-Moore
- Learn from character comparison to skip pointless alignments
- Try alignments in left-to-right order, and try character comparisons in right-to-left order

* Bad character rule
- Upon mismatch skip alignments until a) mismatch becomes a match, or b) P moves past mismatched character

Step1:
T: G C T T C T G C T A C C T T T T G C
P: C C T T T T G C (mismatch at C in T and T in P) (two step adjustment from left to match C)

Step2:
T: G C T T C T G C T A C C T T T T G C
P:       C C T T T T G C (mismatch at A in T and G in P) (This time A does not occur in P) (move all the way past the mismatched character)

Step3:
T: G C T T C T G C T A C C T T T T G C
P:                     C C T T T T G C

* Good suffix rule
- Let t= substring matched by inner loop, skip until a) there are no mismatches between P and t or b) P moves past t

Step1:
               | t |
T: G C T G C C T A C T T A C T T A C T T A C T T A C G C G A A
P: C T T A C T T A C

Step2:
               |     t     |
T: G C T G C C T A C T T A C T T A C T T A C T T A C G C G A A
P:         C T T A C T T A C

Step 3:
T: G C T G C C T A C T T A C T T A C T T A C T T A C G C G A A
P:                 C T T A C T T A C
'''

####################################
# Boyer-Moore Putting all together #
####################################

'''
Whenever we have mismatch we are going to try both rules

Use bad character or good suffix rule, whichever skips more

Step1:
T: G T T A T A G C T G A T C G C G G C G T A G C G G C G A A
P: G T A G C G G C G (no matching character, can't do good suffix rule) (bc:6, gs:0 go with bad character rule)

Step2:
T: G T T A T A G C T G A T C G C G G C G T A G C G G C G A A
P:               G T A G C G G C G (bad character rule does not skip any alignment) (bc:0, gs:2, good suffix)

Step3:
T: G T T A T A G C T G A T C G C G G C G T A G C G G C G A A
P:                     G T A G C G G C G (bc:2, gs:7, good suffix)

Step4:
Step3:
T: G T T A T A G C T G A T C G C G G C G T A G C G G C G A A
P:                                     G T A G C G G C G (match)

Here we skipped 15 alignments, faster than naive matching

Booyer-Moore: Preprocessing
- Pre-calculate skips: For bad chacter rule, P=TCGC
- Because Booyer-Moore reference table for skip numbers
- Only pattern P is needed to build this table, not text T

Table may follow this form
  T C G C
A 0 1 2 3
C 0 - 0 -
G 0 1 - 0
T - 0 1 2

- So we must first build lookup table
'''

##################################
# Diversion: Repetitive elements #
##################################

'''
Real genomes are very unlike "random" genomes
Human genomes is extremely repetitive from transposome (45% of genomes are from transposable elements)

This repeats create ambiguity and problematics in read alignment
'''

############################
# Implementing Boyer-Moore #
############################

'''
Load from jupyter notebook
'''

import string

def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)
    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: Zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]

'''
First Demonstration of Bad character and Good suffix rules
'''
# Bad Character rule
# T = GCTAGCTC
# P = TCAA
p='TCAA'
p_bm=BoyerMoore(p)
p_bm.bad_character_rule(2,'T')

# Good suffix rule
# T = GCTAGCTC
# P = ACTA

p='ACTA'
p_bm=BoyerMoore(p)
p_bm.good_suffix_rule(0) # as 0 position only mismatches

# Match skip
# ACACGCTC
# ACAC
p= 'ACAC'
p_bm=BoyerMoore(p)
p_bm.match_skip()

'''
Let's implement Boyer-Moore
'''
def boyer_moore(p,p_bm,t):
    i=0
    occurrences=[]
    while i < len(t) - len(p) + 1: # all position without running past of T
        shift=1
        mismatched=False
        for j in range(len(p)-1, -1, -1):
            if not p[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j,t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched=True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs=p_bm.match_skip()
            shift = max(shift, skip_gs)
        i+=shift
    return occurrences

# Let test this
t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'
p_bm=BoyerMoore(p)
boyer_moore(p,p_bm,t)

#################
# Preprocessing #
#################

'''
Boyer-Moore build pre-process the pattern P in order to build lookup table that help it use the bc and gs rules

Naive algorithm did not use pre-processing
Input: P and T -> Results

With Boyler-Moore
Make lookup tables for bad characters & good suffix rules

In matching problem if text T is same with varying P, it make sense to make lookup tables for T

Algorithm that preprocesses T is offline, otherwise is online

online: Naive algorithm, Boyer-Moore
offline: web search engine, read alignment
'''

################################
# Indexing and the k-mer index #
################################

'''
Offline Algorithm, a method that pre-process T

Index
Analogy1
- As an analogy when you read book, one way to pre-process it is by indexing it.
- Key terms ordered alphabetically, with associated page #s
Analogy2
- Grocery store items grouped into aisles

Indexing DNA

T: C G T G C G T G C T T

Index of T:
C G T G C : 0, 4
G C G T G : 3
G T G C G : 1
G T G C T : 5
T G C C T : 2
and so on

Here, we used k mer, substring of length k. For example, above is 5-mer index

If P: G C G T G C
How we query index?
- First five character in P is seen at index 3 of T

After this one more step since we are not fully matched. This additional work is called "verification".

What if we take second 5 letter from P? None of them will fail to match
But here we have match at 0, 4 so we need to do verification.

If n-mer matches, we call it "index hit", but not all index hit lead to matches
'''

##################################
# Ordered structure for indexing #
##################################

'''
We discussed how to build and query a k-mer index, an index that's built by taking all the k-mers of the text T
and adding them to a data structure that maps each k-mer to a list of all the offsets where it occurred in the text
is called "multimap"

A k-mer could occur many places within the genome

1. database based on ordering
- we first get list of 3-mer and order them alphabetically
- To query, we do binary search
- Do bisection repeatedly until match
- log2(n) bisections per query

Python provide binary search module
- bisect.bisect_left(a,x): Leftmost offset where x can be inserted into a to maintain order
'''

###########################
# Hash Table for indexing #
###########################

'''
2. Hash table
- Consist of array of buckets
- Hash function h maps 3-mers to buckets
- Bucket entry consists of key: 3-mer, value: offset
- Some distinct 3-mer end up in same bucket, and called "collision"

To query this hash table, we look for bucket that contains 3-mer in Hash table

Python dictionary type is implemtation of hash table
'''

##############################
# Implementing a k-mer index #
##############################

import bisect

class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t)-k+1):
            self.index.append((t[i:i+k],i))
        self.index.sort()

    def query(self, p):
        kmer = p[:self.k]
        i=bisect.bisect_left(self.index, (kmer,-1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i+=1
        return hits

def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k: i+len(p)]:
            offsets.append(i) # verification
    return offsets

t= 'GCTACGATCTAGAATCTA'
p= 'TCTA'

index=Index(t,2)
print(queryIndex(p,t,index))

############################
# Variation on k-mer index #
############################

'''
Two variation on the theme of k-mer indexes

What if we take k-mer at even offset?
- Advantage is that the index is smaller, faster to query
- But since we reduced number of index, we need to use first and second kMers

We can build index over sub-sequence
- string of character also occuring in S in the same order
- a given index hit is more likely to lead to a full match

Substrings are also subsequences, subsequences are not necessarily substrings

We can index T in subsequence way, but need to maintain same pattern in drawing, even for P
'''

###################################
# Genome indexes used in research #
###################################

'''
idea about kmer index is simple, we extract every substring of length k from a reference genome
we can also introduce variation

subsequence is also dealt where we extract subsequences from the genome instead of substrings. So subsequences are not necessarily contiguous

Suffix index
- take all the possible suffix and put them in alphabetical order
- to query this, we can use binary search
- but number of index is n(n+1)/2, it can be computationally challenging

To avoid this, we could introduce suffix array of m integers long
1. Suffix tree - used principle of grouping, >= 45GB
2. Suffix array - principle of ordering, >= 12GB
3. FM index - BWT idea is used, ~ 1GB, widely used
'''

###################################################
# Approximate matching, Hamming and edit distance #
###################################################

'''
Approximate matching
- Difference between read and reference occur because of
1. sequencing error
2. natural variation

Due to this algorithm for exact matching is not sufficient, we need algorithm that does approximate matching allowing differences

example of difference
1) mismatch (substitution)
2) insertion
3) deletion

want to be able to describe how different sequences are by "distance"
1. Hamming distance
- For X & Y where |X| = |Y|, minimum # substituions needs to turn one into the other

2. Edit distance
- For X&Y, minimum # edits (substitutions, insertions, and deletions) needed to turn one into the other
- Does not require the sequence in same length
'''

# To apply hamming distance to our naive algorithm
def naiveHamming(p,t,maxDistance):
    occurrences = []
    for i in xrange(len(t)-len(p)+1): # loop over alignment
        nmm = 0
        match=True
        for j in xrange(len(p)):
            if t[i+j] != p[j]:
                nmm+=1
                if nmm>maxDistance:
                    break
        if nmm <= maxDistance:
            occurrences.append(i)
    return occurrences

########################
# Pigeonhole principle #
########################

'''
Wanted: way to apply exact matching algorithm to approximate matching problems

For P, we first divide them into u and v (u and v are two partitions)
with 1 edit
- If P occurs in T with 1 edit, then u and v appears with no edits
with multiple edits
- If P occurs in T with up to k edits, at least one of p1,p2,p3,...,pk+1 must appear with 0 edits

Since you cannot change all k partitions, we can use this principle

Let say we divide p into five partitions, only match we find is p4.
Next thing we have to do is verification step

This algorithm is a bridge between exact matching and approximate matching
'''

#########################################
# Implementing the pigeonhole principle #
#########################################

'''
Approx matching with Boyle-Moore algorithm
'''

def approximate_match(p,t,n):
    segment_length = round(len(p)/(n+1)) # pigeonhole principle divide p by n+1
    all_matches = set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end],alphabet='ACGT')
        matches = boyer_moore(p[start:end],p_bm,t)

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
    return list(all_matches)

p = 'AACTTG'
t = 'CACTTAATTTG'

print(approximate_match(p,t,2))
print(t[5:])
