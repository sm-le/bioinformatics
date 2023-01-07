"""
# String algorithm 7

Given DNA strings in fasta format, return a consensus string and 
composition profile.

Example,
IN: 
    >seq_1
    ATGC
    >seq_2
    AGAC
    >seq_3
    AGGT
OUT: 
    AGGC
    A: 3 0 1 0
    C: 0 0 0 2
    G: 0 2 2 0
    T: 0 1 0 1
"""

import numpy as np
from bioinformatic_stronghold.tools.fasta_reader import FastaFile

def consensus_profile(fasta_file:str) -> str:
    """Extract consensus sequence from
    fasta file and its matrix profile

    Args:
        fasta_file: path to fasta sequence file
    
    Returns:
        str(consensus sequence)
        str(nucleotide score matrix)
    """

    with FastaFile(fasta_file) as fa:
        fasta_entries = fa.collect()

    # Make a list containing 1D array
    fasta_sequences = [np.array(list(i['sequence']), dtype=str) for i in fasta_entries]
    # Make 2D array from the 1D array list
    fasta_sequences = np.stack(fasta_sequences)

    # Iterate all columns to get nucleotide composition
    nucleotide_composition = list()

    for col in range(fasta_sequences.shape[1]):
        sliced_array = fasta_sequences[:,col]
        # Use np.unique to get unique letter and count information
        nuc_comp, nuc_count = np.unique(sliced_array, return_counts=True)
        nucleotide_composition.append({nuc:count for nuc, count in zip(nuc_comp, nuc_count)})

    # Extract consensus sequence and count profile
    consensus_sequence = ""
    nuc_count_status = {"A":"","C":"","G":"","T":""}
    for nucl_comp in nucleotide_composition:
        consensus_sequence += max(nucl_comp, key=lambda key: nucl_comp[key])
        nuc_count_status["A"] += f"{nucl_comp['A']} " if "A" in nucl_comp else f"0 "
        nuc_count_status["C"] += f"{nucl_comp['C']} " if "C" in nucl_comp else f"0 "
        nuc_count_status["G"] += f"{nucl_comp['G']} " if "G" in nucl_comp else f"0 "
        nuc_count_status["T"] += f"{nucl_comp['T']} " if "T" in nucl_comp else f"0 "
    
    return f"{consensus_sequence}\n" \
           f"A: {nuc_count_status['A']}\n" \
           f"C: {nuc_count_status['C']}\n" \
           f"G: {nuc_count_status['G']}\n" \
           f"T: {nuc_count_status['T']}"