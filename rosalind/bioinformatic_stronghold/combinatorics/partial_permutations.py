"""
# combinatorics 4

Given positive integer n and k, return the 
total number of partial permutation modulo 1M.

Example,
IN:
    21 7
OUT: 
    51200
"""

from math import perm

def partial_permutation(n:int, k:int) -> int:
    """Total number of partial permutations of 
    k object can be formed from n objects.
    
    n! = n(n-1) ... (2) (1)
    n! / (n-k) ! = partial permutation 
    """

    pp = perm(n, k)

    return pp % 1000000