"""
# String algorithm 2

A DNA string 't' corresponds to 'u' upon transcription.

Given a DNA string 't' with length ~1000nt, transcribe DNA to RNA.

Example,
IN: GATGGAACTTGACTACGTAAATT
OUT: GAUGGAACUUGACUACGUAAAUU
"""

def transcribe_dna_v1(sequence:str) -> str:
    """Transcribe DNA to RNA (T -> U)

    Args:
        sequence: a nucleotide sequence only comprised of ACGT
    
    Returns:
        "ACGU..."
    """

    sequence = sequence.upper().replace("T","U")

    return sequence

def transcribe_dna_v2(sequence:str) -> str:
    """Transcribe DNA to RNA (T -> U)

    Args:
        sequence: a nucleotide sequence only comprised of ACGT
    
    Returns:
        "ACGU..."
    """

    transcribe_table = {
                            "A":"A",
                            "C":"C",
                            "G":"G",
                            "T":"U"
                    }

    sequence = "".join([transcribe_table[i] for i in sequence])

    return sequence

def transcribe_dna_v3(sequence:str) -> str:
    """Transcribe DNA to RNA (T -> U)

    Args:
        sequence: a nucleotide sequence only comprised of ACGT
    
    Returns:
        "ACGU..."
    """

    transcribe_table = {
                            "A":"A",
                            "C":"C",
                            "G":"G",
                            "T":"U"
                    }

    sequence = "".join(map(lambda x: transcribe_table[x], sequence))

    return sequence