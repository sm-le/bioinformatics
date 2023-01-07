"""
# alignment 2

Given two DNA strings, return the transition/transversion ratio

Example,
IN:
    >fa_1
    GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
    AGTACGGGCATCAACCCAGTT
    >fa_2
    TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
    GGTACGAGTGTTCCTTTGGGT
OUT: 1.21428571429
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

tvt_matrix = {
                "A" : {
                        "A" : "no change",
                        "G" : "transition",
                        "T" : "transversion",
                        "C" : "transversion"

                    },
                "C" : {
                        "A" : "transversion",
                        "G" : "transversion",
                        "T" : "transition",
                        "C" : "no change"

                    },
                "G" : {
                        "A" : "transition",
                        "G" : "no change",
                        "T" : "transversion",
                        "C" : "transversion"

                    },
                "T" : {
                        "A" : "transversion",
                        "G" : "transversion",
                        "T" : "no change",
                        "C" : "transition"

                    }
            }


def transition(i:str, j:str) -> bool:
    """Return True or False on transition

    Args:
        i: subject nucleotide 1
        j: subject nucleotide 2
    
    Returns:
        True or False on transition
    """

    if tvt_matrix[i][j] == "transition":
        return True
    else:
        return False

def transversion(i:str, j:str) -> bool:
    """Return True or False on transversion

    Args:
        i: subject nucleotide 1
        j: subject nucleotide 2
    
    Returns:
        True or False on transversion
    """

    if tvt_matrix[i][j] == "transversion":
        return True
    else:
        return False

def count_transition_transversion(sequence_1:str, sequence_2:str) -> float:
    """Calculate transition transversion ratio between 
    two DNA sequences and return the ratio

    Args:
        sequence_1: first nucleotide sequence
        sequence_2: second nucleotide sequence

    Returns:
        The transition transversion ratio
    """

    ctransition = sum([transition(i,j) for i, j in zip(list(sequence_1), list(sequence_2))])
    ctransversion = sum([transversion(i,j) for i, j in zip(list(sequence_1), list(sequence_2))])

    return ctransition / ctransversion

def transition_transversion(fasta_file:str) -> float:
    """Read fasta and calculate transition transversion 
    ratio between two DNA sequences and return the ratio
    
    Args:
        fasta_file = fasta file
    Returns:
        The transition transversion ratio
    """

    with FastaFile(fasta_file) as fa:
        records = fa.collect()

    tt_ratio = count_transition_transversion(records[0]["sequence"], records[1]["sequence"])

    return tt_ratio
