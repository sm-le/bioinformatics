# 1. The hidden in DNA is similar to gold-bug assumption. We find the most frequent k-mer within a set of words. 

    #- k-mer refer to a string of length k and defined as Count(Text,Pattern)
    #- we should account overlapping occurrence of Pattern in Text. -> use sliding window

# Pattern - most frequent k-mer in Text if maximize Count.

#Count ...

T = "ATTATATTATATTTATTTATTGCAGCTAGGTACGTGCACACGGTTTACCGGTAGGCACGGT"


def PatternCount(Text, Pattern):
    count = 0
    for i in range(0, (len(Text)-len(Pattern))):
        if Text[i:i+len(Pattern)] == Pattern:
            count+=1
    return count
print(PatternCount(T,"AT"))

# A heuristic algorithm to find k-mer is to find the most frequent k-mer in a string called Text.
# The above algorithm finds a certain k-mer occuring within "Text", however, it does not count the most frequent k-mer happeing within the string.
# To solve this problem, we first need to implement "FrequentWords" that generates array count where count(i) stores Count(Text, Pattern) where Pattern = Text(i,k)

def FrequentWords(Text, k):
    FrequentPatterns = list()
    count = list()
    for i in range(0, len(Text)-int(k)):
        pattern = Text[i:i+k]
        count.append(PatternCount(Text,pattern))
    maxcount = max(count)
    for i in range(0, len(Text)-int(k)):
        if count[i] == maxcount:
            FrequentPatterns.append(Text[i:i+k])
    
    # remove duplicates from list
    FrequentPatterns = set(FrequentPatterns)
    
    return FrequentPatterns

print(FrequentWords(T, 2))

# Although this finds the most frequent k-mer in the text, it is not very efficient in large text. This algorithm must call (Text-k+1) * (Text-k+1) * k step to complete execution. 
# Often this is refer to O(|Text|^2,*k) complexity. 

# 2. Some hidden messages are more elusive than others

    # Additional occurrences can be found in sequences which have one nucleotide difference.
    # This allows fitting to slightly varied target.

# When position p_{i} and q_{i} mismatches for sequences p_{1.2.3.....k} and q_{1.2.3.....k}, we call number of mismatches between p and q is Hamming distance.
# or HammingDistance(p,q)

def HammingDistance(p1,p2):
    count = 0
    for i in range(len(p1)):
        if p1[i] != p2[i]:
            count += 1
    return count

def ApproximatePatternCount(Text,Pattern,d):
    count = 0
    for i in range(0, len(Text)-len(Pattern)):
        Psub = Text[i:i+len(Pattern)]

        if HammingDistance(Pattern,Psub) <= 1:
            count += 1
    return count

print(ApproximatePatternCount("AAAATCGTGCGTTAATAACGGTGCGAAAAC","AAAAA",1))

# multiple local minima in skew graph theory

    # 1. multiple OriC
    # 2. genome rearrangement
    # 3. horizontal gene transfer₩₩