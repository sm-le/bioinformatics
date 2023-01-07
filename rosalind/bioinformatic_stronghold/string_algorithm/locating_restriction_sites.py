"""
# String algorithm 9

Given a DNA string in fasta format, return position and 
length of every reverse palindrome which size between 4 and 12.

Example,
IN:
    >fa_1
    TCAATGCATGCGGGTCTATATGCAT
OUT: 
    4 6
    5 4
    6 6
    7 4
    17 4
    18 4
    20 6
    21 4
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

def reverse_complement_dna(dna_sequence:str) -> str:
    """Reverse complement DNA

    Args:
        dna_sequence: target dna sequence

    Returns:
        str(DNA)
    """

    reverse_comp_table = {
                                    "A":"T",
                                    "C":"G",
                                    "G":"C",
                                    "T":"A"
                            }

    dna_sequence = ''.join(map(lambda x: reverse_comp_table[x], dna_sequence[::-1]))

    return dna_sequence

def palindrome_finder(fasta_file:str) -> list:
    """Find palindrome and return position and
    length

    Args:
        fasta_file: target fasta file
    
    Returns:
        list([pos, length])
    """
    
    with FastaFile(fasta_file) as fa:
        records = fa.collect()

    palindrome_summary = dict()

    for record in records:
        seq = record["sequence"]
        rseq = reverse_complement_dna(seq)

        palindrome = list()
        for lnx in range(4, 13):
            for idx in range(len(seq)-lnx+1):
                if seq[idx:idx+lnx] == rseq[len(rseq)-idx-lnx:len(rseq)-idx] and len(seq[idx:idx+lnx]) >= 4:
                    palindrome.append([idx+1, lnx])

        palindrome_summary[record["description"]] = palindrome
    
    return palindrome_summary