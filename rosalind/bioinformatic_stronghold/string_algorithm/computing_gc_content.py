"""
# String algorithm 4

A GC content of a DNA string is given by the percentages of G and C in the string.

Given 10 fasta formatted sequences, return id and gc-content with the highest gc content.

Example,
IN: 
    >fasta_1
    ATATC
    >fasta_2
    GGTCA
    >fasta_3
    GGCCT
OUT:
    fasta_3:80.0
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

def highest_gc_v1(fasta_path:str) -> str:
    """Calculate GC content of each sequence 
    and extract sequence header and GC percentage 
    of one with the highest GC content

    Args:
        fasta_path: path to fasta file

    Returns:
        string({fasta_header}\n{gc_content})
    """

    with FastaFile(fasta_path) as fa:
        fasta_entries = fa.collect()

    highest_gc = sorted(fasta_entries, key= lambda x: x["gc"], reverse=True)[0]
    
    return f"{highest_gc['description']}: " \
           f"{highest_gc['gc'] * 100}"