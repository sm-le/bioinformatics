"""
# proteomics 1

Given a list of UnitProt ID, return location of N-glycosylation motif. 

Example,
IN:
    B5ZC00
OUT:
    85 118 142 306 395
"""

import bioinformatic_stronghold.tools.uniprot_parser as upp
import re

def protein_motif_finder_v1(motif_name:str, accession:str) -> str:
    """Find targeted protein motif and return 
    all occurrences 

    Args:
        motif_name: name of motif to search
        accession: targeted uniprot accession 
    
    Returns:
        str(locations of targeted motif)
    """

    # Read patterns
    if motif_name == "N-gly":
        motif_pattern = r'N[^P][ST][^P]'

    # Uniprot download and parse
    uniprot = upp.parse_uniprot(accession)
    
    
    # dictionary to save result
    motif_location = dict()

    # Find location
    locs = [nuc.start()+1 for nuc in re.finditer(motif_pattern, uniprot["sequence"])]
    # Save result
    motif_location["header"] = accession
    motif_location["location"] = locs

    return motif_location

def protein_motif_finder_v2(motif_name:str, accession:str) -> str:
    """Find targeted protein motif and return 
    all occurrences 

    Args:
        motif_name: name of motif to search
        accession: targeted uniprot accession 
    
    Returns:
        str(locations of targeted motif)
    """

    # Read patterns
    if motif_name == "N-gly":
        motif_pattern = r'[N][^P][ST][^P]'

    # Uniprot download and parse
    uniprot = upp.parse_uniprot(accession)
    
    # dictionary to save result
    motif_location = dict()

    # Find location
    ## silenced find iter due to lack of adjacent location finder
    # locs = [nuc.start()+1 for nuc in re.finditer(motif_pattern, uniprot["sequence"])]
    match_pos = 0
    locs = list()
    while True:
        try:
            match_info = re.search(motif_pattern, uniprot["sequence"][match_pos:]).start()
            match_pos += match_info + 1
            locs.append(match_pos)
        except:
            break

    # Save result
    motif_location["header"] = accession
    motif_location["location"] = locs

    return motif_location