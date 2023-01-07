"""
# combinatorics 1

Given a max:1000aa protein sequence, return total probable number of 
RNA sequence module 1M. 

Example,
IN:
    MA
OUT: 
    12
"""

def possible_mRNA_from_protein(target_peptide):
    """Calculate all possible number of 
    mRNA inferred from a sequence of protein
    
    Args:
        target_peptide: target sequence of amino acid
    
    Returns:
        int(combination modulo a million)
    """

    aa_table = {
                "F":["UUU","UUC"],
                "L":["UUA","UUG","CUU","CUC","CUA","CUG"],
                "I":["AUU","AUC","AUA"],
                "M":["AUG"],
                "V":["GUU","GUC","GUA","GUG"],
                "S":["UCU","UCC","UCA","UCG","AGU","AGC"],
                "P":["CCU","CCC","CCA","CCG"],
                "T":["ACU","ACC","ACA","ACG"],
                "A":["GCU","GCC","GCA","GCG"],
                "Y":["UAU","UAC"],
                "*":["UAA","UAG","UGA"],
                "H":["CAU","CAC"],
                "Q":["CAA","CAG"],
                "N":["AAU","AAC"],
                "K":["AAA","AAG"],
                "D":["GAU","GAC"],
                "E":["GAA","GAG"],
                "C":["UGU","UGC"],
                "W":["UGG"],
                "R":["CGU","CGC","CGA","CGG","AGA","AGG"],
                "G":["GGU","GGC","GGA","GGG"]
                }

    target_peptide += "*"

    combinations = 1
    for aa in target_peptide:
        combinations = combinations * len(aa_table[aa])

    return combinations % 1000000