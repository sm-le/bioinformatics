"""
# Probability 4

Given a DNA string and an array, 
return an array having the same length as given 
array where element of the array represent common 
logarithm of the probability

Example,
IN: ACGATACAA
    0.129 0.287 0.423 0.476 0.641 0.742 0.783
OUT: 
    -5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009
"""

import math

def random_string_probability(dna:str, array:str) -> str:
    """Compute probability of kth random dna string from 
    array[k] match to dna

    Args:
        dna: a dna string
        array: GC probability array
    """

    array = list(map(float, array.split(" ")))

    rs_probability = list()

    for pGC in array:

        nP = {
                "G": pGC / 2,
                "C": pGC / 2,
                "A": (1-pGC) / 2,
                "T": (1-pGC) / 2
            }

        em_prob = list(map(lambda x: nP[x], dna))
        em_log = round(math.log(math.prod(em_prob),10),3)

        rs_probability.append(em_log)

    return " ".join(map(str,rs_probability))