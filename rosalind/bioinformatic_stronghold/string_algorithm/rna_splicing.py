"""
# String algorithm 10

Given a DNA string and a collection of substring in fasta format, 
return a protein string from an exon.

Example,
IN: 
    >fasta_1
    ATGGAGTTC
    >fasta_2
    GAG
OUT: 
    MF
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

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

def translate_rna(rna_sequence:str) -> list:
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

    aa_sequence = ''.join(map(lambda x: codon_table[x], [rna_sequence[i:i+3] for i in range(0,len(rna_sequence),3)]))

    return aa_sequence

def translate_exon(fasta_file:str):
    """Remove intron from dna and translate
    
    Args:
        fasta_file: target fasta
    
    Returns:
        str(aa)
    """

    with FastaFile(fasta_file) as fa:
        records = fa.collect()

    dna = records[0]["sequence"]
    introns = records[1:]

    for intron in introns:
        dna = "".join(dna.split(intron["sequence"]))

    aa = translate_rna(transcribe_dna(dna))

    return aa