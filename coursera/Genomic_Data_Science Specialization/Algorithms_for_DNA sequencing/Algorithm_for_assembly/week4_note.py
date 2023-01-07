# @sungml92
# Algorithm for DNA sequencing
# week 4: Algorithm for assembly

################
# Introduction #
################

'''
First discussion will be on assembly problem called the shortest common super strength problem

Second we will look at DeBruin graph to solve when the genome is repetitive

Finally we will discuss some aspect of how real assembly program may work
'''

########################################################
# Lecture: The shortest and common superstring problem #
########################################################

'''
Shortest common superstring

- Given set of string S, find SCS(S): shortest string containing the string in S as substrings

e.g) S: BAA AAB BBA ABA ABB BBB AAA BAB
concat(S): BAAAABBBAABAABBBBBAAABAB
SCS(S): AAABBBABAA (we cannot just concat string since we are required to find shortest string)

- This is important when we look at all 6-mer network that we produced in week3
- When we get shorest common string of all 6 6-mers, we will get the synthetic genome
- For this reason, it is important to get genome from the short overlapping string

This give most parsimonous explanation on the reads.

SCS has some downsides

1) NP-complete: no efficient algorithm for large inputs -> not so fast

- Idea, pick order for substrings in S and construct superstring
order 1: AAA AAB ABA ABB BAA BAB BBA BBB
        AAABABBAABABBABBB <- supoerstring 1, if we pick right ordering this will be the shortest super string

and continue doing this, until you try every permutations. So this is very computationally expensive when it comes to many inputs

If S contains n strings, n! (n factorial) orderings are possible. So we called this intractable

'''

#######################################################
# Practical: Implementing shortest common superstring #
#######################################################

'''
Implemeting SCS with brute force
'''

def overlap(a,b,min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start+=1

import itertools

def scs(ss):
    shortest_sup = None
    for ssperm in itertools.permutations(ss): # return every possible reads
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i],ssperm[i+1], min_length=1)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
    return shortest_sup

scs(["ACGGTACGAGC","GAGCTTCGGA","GACACGG"])

###############################################
# Lecture: Greedy shortest common superstring #
###############################################

'''
Brute force is no efficient method of solving problem so it is extremely slow.

Greedy shortest common superstring
- algorithm will make series of decision reducing the length of superstring the most
- we can visualize the greedy shortest common superstring algorithm using an overlap graph
- at each round, we pick longest remaining overlap in the graph
- and merge the node on either side of the edge

Keep in mind that this may result in local minima
'''

##############################################################
# Practical: Implementing greedy shortest common superstring #
##############################################################

'''
first find two reads with maximum overlap and combine them <- repeat this
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

greedy_scs(["ABC","BCA","CAB"],2)

'''
There is discrepancy between greedy_scs and scs as greedy_scs does not always lead to correct shortest superstring
'''

###################################################
# Lecture: Third law of assembly: repeats are bad #
###################################################

'''
In the nature of repetitive genome, finding a shortest superstring for the genome may not incorporate repetitive elements

Third law of assembly
- repeats make assembly difficult

No matter slow or greedy algorithm, they will always collapse repeative copies into fewer copies.

In reality, half the genome is covered by repetitive elements, so this makes assembly very difficult

We will talk about collapse problem
'''

################################################
# Lecture: De Bruijn graphs and Eulerian walks #
################################################

'''
issue with the shortest common superstring formulation, which is that if the genome is repetitive,
it will tend to over collapse the repeats. It will assemble fewer copies of the repeats than are actually there.

De Brujin graph
- it is a directed graph

for example, "tomorrow and tomorrow and tomorrow"
two distinct word is "tomorrow", "and"
every time it directly to certain node new edge will form, no matter how many edges were there previously

In the de brujin graph, after generating kmer of genome it will subsequently generate L/R (k-1)mers
New node will only form if it is unique, for example, AAA (3-mer) will generate AA (L) and AA (R) 2 mers and AA node will generate will self feedback edge

- one edge per k-mer
- one node per distinct k-1-mer

Walking crossing each edge exactly once gives a reconstruction of the genome. This is called Eulerian walk
'''

#########################################
# Practical: Building a De Bruijn graph #
#########################################

def de_brujin_ize(st,k):
    edges = []
    nodes = set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes,edges
nodes, edges = de_brujin_ize("ACGCGTCG",3)
print(nodes)
print(edges)

def visualize_de_brujin(st,k):
    nodes,edges = de_brujin_ize(st,k)
    dot_str = 'digraph "DeBrujin graph" {\n'
    for node in nodes:
        dot_str += ' %s [label="%s"] ;\n' % (node,node)
    for src, dst in edges:
        dot_str += ' %s -> %s ; \n' % (src,dst)
    return dot_str + '}\n'

'''
for use in ipython
%load_ext gvmagic
%dotstr visualize_de_brujin("ACGCGTCG",3)
'''

#########################################
# Lecture: When Eulerian walks go wrong #
#########################################

'''
DeBrujin graph with Eulerian path may solve repetitive genome problem, but itself can encounter many path problem

For example, if we construct a graph on ZABCDABEFABY, k=3 can have two pathways for walking along DeBrujin graph

decreasing kmer length, we are increasing the chance of being affected by repeats since the smaller k means there's more likely to be multiple
occurrences of any given k-mer in the genome which is going to increase the chance that we have multiple different Eulerian walk for the same graph

de Bruijn graph representation is actually a very common and useful way to represent assemblies,
to represent the assembly problem, and in fact, many modern software tools
for assembly use de Bruijn graphs as the internal representation of the relationships between the sequencing reads.
So while shortest common superstring and Eulerian walk are both flawed formulations of the assembly problem,
the overlap graph and the de Bruijn graph are both going to continue to be very useful to us in practice.
'''

###################################
# Lecture: Assemblers in practice #
###################################

'''
Two kind of graphs
1. Overlap graph -> Overlap-Layout-Consensus (OLC) assembly
2. De Brujin graph -> De Brujin Graph based (DBG) assembly
3. Common thing is that we are not using SCS or Eulerian walk

First step is to build a graph, but "graph will be very messy"
Reason of messy graph
1. Sequencing error leads to dead ends. True reconstruction is to ignore these dead ends
2. Transitively inferable edges
3. Polyploidy
    - to avoid this, collapse the graph and leave a note
4. Repeats

We deal with this by chopping assembly into pieces, so we can put some portion of puzzle together
and partial reconstruction that we can put together are called contigs

So after product is a set of contigs, not a single assembled genome
'''

################################
# Lecture: The future is long? #
################################

'''
The longer the reads are the more likely we are to get a read that anchors
some repetitive sequence that glues it with some surrounding non-repetitive sequence.
And that's what tells us where the repetitive sequence should go in the assembly.

1. paired-end sequencing

2. single molecule sequencer
'''

##############################################
# Lecture: Computer science and life science #
##############################################
