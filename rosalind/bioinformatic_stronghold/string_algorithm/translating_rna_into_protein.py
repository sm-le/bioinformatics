"""
# String algorithm 5

Given an RNA string 's', return translated protein string.

Example,
IN: AUGGCC
OUT: MA
"""

def translate_rna(sequence:str) -> str:
    """Translate RNA to protein

    Args:
        sequence: a nucleotide sequence only comprised of ACGU...
    
    Returns:
        "amino acid string"
    """

    codon_table = {
                    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
                    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
                    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
                    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
                    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
                    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
                    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
                    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
                    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
                    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
                    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
                    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
                    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
                    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
                    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
                    "GGU":"G","GGC":"G","GGA":"G","GGG":"G"
                }

    sequence = "".join([codon_table[sequence[i:i+3]] for i in range(0,len(sequence),3)])

    return sequence