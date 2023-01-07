"""
# String algorithm 6

Given two DNA strings 's' and 't', find all locations of 't' in 's'

Example,
IN: 
    GAAGAA
    GAA
OUT: 
    1 4
"""

def motif_finder_v1(query:str, subject:str) -> str:
    """Find subject motif within query

    Args:
        query: query string
        subject: a subject to find within a query string
    
    Returns:
        str(start locations)
    """

    return " ".join([str(i+1) for i in range(len(query)) 
                    if query[i:i+len(subject)] == subject])

def motif_finder_v2(query:str, subject:str) -> str:
    """Find subject motif within query

    Args:
        query: query string
        subject: a subject to find within a query string
    
    Returns:
        str(start locations)
    """

    return " ".join(list(set([str(query[i:].find(subject)+1+i) 
                    for i in range(0, len(query)-len(subject)) 
                    if query[i:].find(subject) != -1])))