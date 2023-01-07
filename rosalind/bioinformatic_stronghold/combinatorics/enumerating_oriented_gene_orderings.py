"""
# combinatorics 5

Given a positive integer, return total number of 
signed permutations

Example,
IN:
    2
OUT: 
    8
    -1 -2
    -1 2
    1 -2
    1 2
    -2 -1
    -2 1
    2 -1
    2 1
"""

import itertools

l = 5

def signed_permutation(length:int) -> list:
    """Make permutation of given length 
    including both sign + and -

    Args:
        length: length of permutation
    Returns:
        list(permutation of length: length)
    """

    # make permutations
    perms = list(map(list,itertools.permutations(list(range(1,length+1)))))

    # make empty list to save signed permutation
    signed_perms = list()

    # iterate each permutation and multiply signed array
    for perm in perms:
        for signed_array in itertools.product([-1,1],repeat=len(perm)):
            signed_perms.append([p*array for p,array in zip(perm,signed_array)])

    return signed_perms