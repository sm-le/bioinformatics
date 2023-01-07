"""
# introduction 5

Given a string with letters, return a corpus

Example,
IN:
    yes you can YOU CAN
OUT: 
    yes 1
    you 1
    can 1
    YOU 1
    CAN 1
"""

def string_to_dictionary_v1(string:str) -> str:
    """Read a string and make a word
    corpus

    Args:
        string: a target string

    Returns:
        a count dictionary of every letter
    """
    
    my_corpus = dict()

    for letter in string.split(" "):
        if letter in my_corpus:
            my_corpus[letter] += 1
        else:
            my_corpus[letter] = 1
    
    return my_corpus
