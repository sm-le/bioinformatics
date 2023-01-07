"""
# Probability 1

Given a population k + m + n comprising three positive integer for homozygous 
dominant (k), heterozygous (m), and homozygous recessive (n) comprising, return
probability that two randomly selected mating organism will produce a dominant 
allele.

Example,
IN: 2 2 2
OUT: 0.78333
"""

def dominant_allele(k:int, m:int, n:int) -> float:
    """Calculcate probability of dominant allele on 
    random mating within a sample population
    
    Args:
        k: a number of sample with homozygous dominant
        m: a number of sample with heterzygous
        n: a number of sample with homozygous
    
    returns:
        a probability of a offspring with dominant allele
    """

    """
    Pr(dominant allele) = 1 - Pr(recessive allele)
    Pr(recessive allele) = Pr(X=n and X=n) + Pr(X=m and X=n) * (1/2) +
                           Pr(X=n and X=m) * (1/2) + Pr(X=m and X=m) * (1/4)

    Pr(X=n and X=n) = (( n / ( k+m+n )) * (( n-1 ) / ( k+m+n-1 )))
    Pr(X=m and X=n) = (( n / ( k+m+n )) * (( m ) / ( k+m+n-1 ))
    Pr(X=n and X=m) = (( m / ( k+m+n )) * (( n ) / ( k+m+n-1 ))
    Pr(X=m and X=m) = (( m / ( k+m+n )) * (( m-1 ) / ( k+m+n-1 ))

    """

    return 1 - ((( n / ( k+m+n )) * (( n-1 ) / ( k+m+n-1 ))) + \
           (( n / ( k+m+n )) * (( m ) / ( k+m+n-1 )) * ( 1/2 )) + \
           (( m / ( k+m+n )) * (( n ) / ( k+m+n-1 )) * ( 1/2 )) + \
           (( m / ( k+m+n )) * (( m-1 ) / ( k+m+n-1 )) * ( 1/4 )))