"""
# String algorithm 1

Given a DNA 's' of length ~ 1000nt, count 'A','C','G','T' within s.

Example,
IN: AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
OUT: 20 12 17 21
"""

def count_ACGT_v1(sequence:str) -> str:
    """Count ACGT individually given input 
    nucleotide sequence

    Args:
        sequence: a nucleotide sequence

    Returns:
        string(count(A)\scount(C)\scount(G)\scount(T))
    """

    # First direct counting
    A = sequence.count("A")
    C = sequence.count("C")
    G = sequence.count("G")
    T = sequence.count("T")

    return f"{A} {C} {G} {T}"

def count_ACGT_v2(sequence:str) -> str:
    """Count ACGT individually given input 
    nucleotide sequence

    Args:
        sequence: a nucleotide sequence

    Returns:
        string(count(A)\scount(C)\scount(G)\scount(T))
    """

    # Loop dictionary
    nuc_count = dict()
    for nuc in sequence:
        try:
            nuc_count[nuc] += 1
        except:
            nuc_count[nuc] = 1

    return f"{nuc_count['A']} {nuc_count['C']} {nuc_count['G']} {nuc_count['T']}"

def count_ACGT_v3(sequence:str) -> str:
    """Count ACGT individually given input 
    nucleotide sequence

    Args:
        sequence: a nucleotide sequence
    
    Returns:
        string(count(A)\scount(C)\scount(G)\scount(T))
    """
    from collections import Counter

    # Use counter
    sequence = list(sequence)
    nuc_count = Counter(sequence)

    return f"{nuc_count['A']} {nuc_count['C']} {nuc_count['G']} {nuc_count['T']}"