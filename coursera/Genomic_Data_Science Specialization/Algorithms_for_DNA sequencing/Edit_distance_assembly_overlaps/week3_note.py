# @sungml92
# Algorithm for DNA sequencing
# Week 3: edit distance, assembly, overlaps

################
# Introduction #
################

'''
1. discussing dynamic programming for edit distance, and variation on the algorithm
2. genome assembly from a scratch
'''

##############################################
# Lecture: Solving the edit distance problem #
##############################################

'''
Dynamic programming algorithm does not depend on exact matching algorithm, so it is a separate kind.

Two different method in measuring distance
    1. Hamming distance is the distance where two strings were equal length and number of substitution required to turn one to the other
    2. Edit distance is the distance where minimum number of edits (substitution, insertions, deletions) needed to turn one to the other
'''
# It is easy to find hamming distance
def hammingDistance(x,y):
    nmm = 0
    for i in xrange(0,len(x)):
        if x[i] != y[i]:
            nmm+=1
    return nmm

# But edit distance is not that easy to construct
'''
Relationship between hamming distance and edit distance
    - editdistance <= hammigdistance
    - ability to use insertion and deletion lead us to get fewer changes than hamming distance

If x and y are different length, what can we say about editdistance(x,y)?
    - editDistance(X,Y) >= ||X|-|Y||
    - at very least we need to make the same length

If we know editdistance between two prefixes it is easier to compute editdistance
    - editDistance(X[:-1],Y[:-1])

Let's called prefixes of each string as alpha and beta
so that strings can be represented as
    alpha x
    beta y
and therefore, edist(ax,by) = min of edist(a,b)+delta(x,y), edist(ax,b)+1, edist(a,by)+1
    delta(x,y) = 0 if x==y, else 0 if x!=y

This edit distance can be explained by recursion
'''
# some example can be
def eDistRecursive(a,b):
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    #return len(b) if len(a)==0 else len(a) if len(b)==0
    delt = 1 if a[-1] != b[-1] else 0
    return min(eDistRecursive(a[:-1],b[:-1]) + delt, eDistRecursive(a,b[:-1]) + 1, eDistRecursive(a[:-1],b)+1) # which this exactly correspond to edist(ax,by)

########################################################
# Lecture: Using dynamic programming for edit distance #
########################################################

'''
recursion is very very slow because redundant work is prevalent with recursive coding

Let's simplify the workload by matrix, we will fill the editdist operation in different pair of argument
So that for any pair of prefixes from X & Y, edit distance is calculated once
'''

################################################################
# Practice: Implementing dynamic programming for edit distance #
################################################################

# Dynamic programming
def editDistance(a,b):
    mat = []
    for i in range(len(a)+1):
        mat.append([0]*(len(b)+1)) # make matrix x+1 and y+1 array of 0

    for i in range(len(a)):
        mat[i][0] = i
    for i in range(len(b)):
        mat[0][i] = i

    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            distHor = mat[i][j-1] + 1
            distVer = mat[i-1][j] + 1
            distDiag = mat[i-1][j-1] if a[i-1] == b[j-1] else mat[i-1][j-1]+1

            mat[i][j] = min(distHor, distVer, distDiag)
    return mat[-1][-1]

editDistance("shakr spear", "shakespear")

###################################################
# Lecture: A new solution to approximate matching #
###################################################

'''
Now we want to apply idea of dynamic programming to approximate matching

Imagine the matrix of two alignment we used in editDistance.
 The idea remains the same but two alignment becomes P and T.
  But here we initialize first row with all zeros and column rule does not change. This is because we do not know where in T occur P

To find alignmnet, we must do "tracing back" along the matrix. Re-trace to know how we got to the bottom.

In summary, We fill in a matrix just like this,
and we can look in the final row to detect occurrences of P in T with some number of edits,
and then we can use the trace back to identify the exact location where that match occurs and what the shape of the alignment is.

In practice, edit distance-based approximate matching tends to be used in combination with other techniques
like techniques that we've learned about using an index, using the pigeon hole principle
'''

########################################################
# Lecture: Meet the family: global and local alignment #
########################################################

'''
Lesson on couple variation on the dynamic programming in edit distance and use it to solve global and local alignment problem.

Global alignment
- First we need to understand not all substitutions are equally occuring. We need a penalty system for transition and transversion.
    In reality, transitions are more frequent than transversion.

- Second, indels are less frequent than substitutions.

- For this reason, we need to introduce penalty matrix. This matrix has an element for every kind of penalty that might be incurred in an approximate match.

- for edist logic all we need to change for adding penalty matrices are just adding different contributions.
    Hence, galign(ax,by) = min(galign(a,b)+p(x,y), galign(ax,b)+p(x,-), galign(a,by)+p(-,y))

Local alignment
- Find the most similar pair of substrings from X and Y
- algorithm is similar to global alignment
    Hence, lalign(ax,by) = max(lalign(a,b)+s(x,y), lalign(ax,b)+s(x,-), lalign(a,by)+s(-,y),0)
- Instead, here we use scoring matrix that gives positive score when match, negative score when there are differences.
- With all background of zeros, the values that are non-0 are mostly indicate local alignment region
'''

############################################
# Practical: Implementing global alignment #
############################################

# We are going to modify edit distance to do global alignment

def globalAlignment(a,b):
    alphabet = ['A','C','G','T']
    score = [[0, 4, 2, 4, 8], \
             [4, 0, 4, 2, 8], \
             [2, 4, 0, 4, 8], \
             [4, 2, 4, 0, 8], \
             [8, 8, 8, 8, 8]]
    mat = []
    for i in range(len(a)+1):
        mat.append([0]*(len(b)+1)) # make matrix x+1 and y+1 array of 0

    for i in range(1, len(a)+1):
        mat[i][0] = mat[i-1][0] + score[alphabet.index(a[i-1])][-1] # alphabet.index(x[i-1])] what row to look at our score in the array and take last value
    for i in range(1, len(b)+1):
        mat[0][i] = mat[0][i-1] + score[-1][alphabet.index(b[i-1])]

    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            distHor = mat[i][j-1] + score[-1][alphabet.index(b[j-1])]
            distVer = mat[i-1][j] + score[alphabet.index(a[i-1])][-1]
            distDiag = mat[i-1][j-1] if a[i-1] == b[j-1] else mat[i-1][j-1]+score[alphabet.index(a[i-1])][alphabet.index(b[j-1])]

            mat[i][j] = min(distHor, distVer, distDiag)
    return mat[-1][-1], mat

# Let's test this
a='TCTCTACTCATC'
b='TCTCTACTCATC'
globalAlignment(a,b)
b='TCTCTAAACATC'
globalAlignment(a,b)

########################################
# Lecture: Read alignment in the field #
########################################

'''
Index and dynamic programming work really well together

Index allows us to rapidly home a small set of candidate location in the genome where there is a "hit", may as well act as a filter.
    but indexing is not great when it comes to approximate matching

Instead, we can directly do dynamic programming and skip indexing. But matrix become really big.
That is why we need index + dynamic programming.
'''

###########################################
# Lecture: Assembly: working from scratch #
###########################################

'''
read alignment is where we are putting puzzle in the complete picture we know each piece belong to

In de novo shotgun assembly, we do not have the benefit of complete picture.
- Fundamentally difficult and computationally challenging.

There are confounding tools to help tackling this problem
'''

##############################################
# Lecture: First and second laws of assembly #
##############################################

'''
Assembly is recontruction of a sequence from multiple reads

Terms in assembly

1. Coverage
    - Coverage is amount of redundant information we have about the genome
    - Used as an evidence of base occuring at the position
    - But not always has same composition. e.g. mixture of a and g
    - We can also calculate overall coverage


Law of assembly

1. First law of assembly
- If a suffix of read A is similar to a prefix of read B,
    then A and B might overlap in the genome. So we can glue them together

- Sometimes there is a mismatch between overlapping read
    1. sequencing errors
    2. polyploidy, human has 2 copies of each chromosome, and copies can differ

2. Second law of assembly
- More coverage leads to more and longer overlaps
'''

###########################
# Lecture: Overlap graphs #
###########################

'''
overlap is a glue that assemble genome back together

Here we will learn to represent overlaps altogether.
 For this, we will draw directed graph by nodes and edges

We can make directed graph to represent all overlap relationship
- where each node is read and draw edge A->B when suffix of A overlaps prefix of B
- connect node with overlaps greater than cut-off to retrieve full genome
'''

##############################################
# Practical: Overlaps between pairs of reads #
##############################################

# function to find overlaps between two strings

def overlap(a,b,min_length=3):
    start=0

    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start+=1
overlap("TTAGCT","TAGCTA",min_length=3)

####################################################
# Practical: Finding and representing all overlaps #
####################################################

# extend overlap function to create naive overlap map
from itertools import permutation

def naive_overlap_map(reads,k):
    olaps = {}
    for a,b in permutations(reads,2): # all pairs of read
        olen = overlap(a,b,min_length=k)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps
