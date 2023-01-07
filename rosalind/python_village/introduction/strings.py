"""
# introduction 2

Given a text and four positive integers a, b, c, and d, return
sliced string of a to b and c to d.

Example,
IN: hahahahahaha
    0 2 2 4
OUT: ha ha
"""

def string_slicer(string:str,a:int,b:int,c:int,d:int) -> str:
    """Slice string given a, b, c, d

    Args:
        string: a target string to slice
        a: first position of a first
        b: second position of a first
        c: first position of a second
        d: second position of a second

    Returns:
        two sliced strings
    """

    return string[a:b+1], string[c:d+1]