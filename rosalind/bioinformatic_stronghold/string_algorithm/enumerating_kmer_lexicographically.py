"""
# String algorithm 11

Given a DNA symbol and integer, return all string with defined length

Example,
IN: 
    A C G T
    2
OUT: 
    AA
    AC
    AG
    AT
    CA
    CC
    CG
    CT
    GA
    GC
    GG
    GT
    TA
    TC
    TG
    TT
"""

from itertools import product

def k_product(letter_set:str, k:int):
    """Generate product of k length from a 
    letter set

    Args:
        letter_set: nucleotide letter set
    
    Returns:
        list(k length permutations)
    """

    k_prod = sorted(list(product(list(letter_set), repeat=k)))

    return k_prod
