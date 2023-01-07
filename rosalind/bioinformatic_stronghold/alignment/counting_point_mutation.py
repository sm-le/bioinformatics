"""
# alignment 1

Given two strings s, t with equal length, a hamming distance is a distance 
between s and t.

Given two DNA strings, return the hamming distance.

Example,
IN:
    TACTAA
    TAATAC
OUT: 2
"""

def hamming_distance_v1(sequence_1:str, sequence_2:str) -> int:
    """Calculate hamming distance between 
    two DNA sequences and return hamming distance 
    aka difference between two string in 
    nucleotide level

    Args:
        sequence_1: first nucleotide sequence
        sequence_2: second nucleotide sequence

    Returns:
        The hamming distance in numeric form
    """

    hamming_distance = sum([i!=j for i, j in zip(list(sequence_1), list(sequence_2))])

    return hamming_distance