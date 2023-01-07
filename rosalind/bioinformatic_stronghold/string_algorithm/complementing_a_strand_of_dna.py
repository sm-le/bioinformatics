"""
# String algorithm 3

A DNA string 'a' complements to 't' so as c and g.

Given a DNA string 's' with length ~1000nt, reverse complement s to s^{c}.

Example,
IN: AAAACCCGGT
OUT: ACCGGGTTTT
"""

def rcomplement_dna_v1(sequence:str) -> str:
    """Complement given dna sequence to reverse
    strand

    Args:
        sequence: a nucleotide sequence

    Returns:
        string("TGCA....")
    """

    complement_table = {
                        "A":"T",
                        "T":"A",
                        "C":"G",
                        "G":"C"
                    }
    
    sequence = "".join(map(lambda x: complement_table[x], sequence[::-1]))

    return sequence

def rcomplement_dna_v2(sequence:str) -> str:
    """Complement given dna sequence to reverse
    strand

    Args:
        sequence: a nucleotide sequence

    Returns:
        string("TGCA....")
    """

    complement_table = {
                        "A":"T",
                        "T":"A",
                        "C":"G",
                        "G":"C"
                    }
    
    sequence = "".join([complement_table[x] for x in  sequence[::-1]])

    return sequence