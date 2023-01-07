"""
# combinatorics 3

Given a positive integer, return total number of 
permutations followed by a list of permutations

Example,
IN:
    3
OUT: 
    6
    1 2 3
    1 3 2
    2 1 3
    2 3 1
    3 1 2
    3 2 1
"""

from itertools import permutations

def n_permutation(number:int) -> list:
    """Generate permutation with 
    given number
    
    Args:
        number: number of element

    Returns:
        list(permutations)
    """

    return list(map(list,permutations(list(range(1,number+1)),number)))