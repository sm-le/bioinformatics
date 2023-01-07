"""
# combinatorics 2

Given a max:1000aa nucleotide sequence, return all candidate protein from 
ORFs.

Example,
IN:
    >fa_1
    AATGATG
OUT: 
    MM
    M
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile
import re

def transcribe_dna(dna_sequence:str) -> str:
    """Transcribe DNA to RNA

    Args:
        dna_sequence: target dna sequence

    Returns:
        str(RNA)
    """

    transcribe_table = {
                            "A":"A",
                            "C":"C",
                            "G":"G",
                            "T":"U"
                    }

    rna_sequence = ''.join(map(lambda x: transcribe_table[x], dna_sequence))

    return rna_sequence

def reverse_transcribe_dna(dna_sequence:str) -> str:
    """Transcribe DNA to RNA

    Args:
        dna_sequence: target dna sequence

    Returns:
        str(RNA)
    """

    reverse_transcribe_table = {
                                    "A":"U",
                                    "C":"G",
                                    "G":"C",
                                    "T":"A"
                            }

    rna_sequence = ''.join(map(lambda x: reverse_transcribe_table[x], dna_sequence[::-1]))

    return rna_sequence

def translate_ORF(rna_sequence:str) -> list:
    """Translate RNA to amino acids

    Args:
        rna_sequence: target rna sequence

    Returns:
        str(aa)
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

    
    orf_aa = list()
    read_frame = 3

    for idx in range(read_frame):
        orf = False
        aa = ""

        for i in range(idx, len(rna_sequence)-read_frame-idx,3):
            if codon_table[rna_sequence[i:i+3]] == "M":
                orf = True
            if orf and codon_table[rna_sequence[i:i+3]] == "*":
                orf_aa.append(aa)
                aa = ""
                orf = False
            if orf:
                aa += codon_table[rna_sequence[i:i+3]]
    
    return orf_aa

def retrieve_variations(orf_list:list) -> list:
    """Retrieve all variation within a single orf 

    Args:
        a list of orfs found
    
    Returns:
        list(orfs)
    """

    orf_variations = list()

    for orf in orf_list:
        if orf.count("M") > 1:
            match_pos = 0
            while True:
                try:
                    match_info = re.search("M", orf[match_pos:]).start()
                    match_pos += match_info + 1
                    orf_variations.append(orf[match_pos-1:])
                except:
                    break
        else:
            orf_variations.append(orf)

    return orf_variations

def ORFs(fasta_file) -> dict:
    """Retrieve all opening reading frames 
    from given fasta file

    Args:
        fasta_file: target fasta file
    
    Returns:
        dict(str(fasta_header):dict("fwd":list, "rev":list))
    """

    with FastaFile(fasta_file) as fa:
        records = fa.collect()

    orf_dict = dict()

    for record in records:
        fwd_orfs = retrieve_variations(translate_ORF(transcribe_dna(record["sequence"])))
        rev_orfs = retrieve_variations(translate_ORF(reverse_transcribe_dna(record["sequence"])))
        
        orf_dict[record["description"]] = dict()
        orf_dict[record["description"]]["fwd"] = fwd_orfs
        orf_dict[record["description"]]["rev"] = rev_orfs

    return orf_dict