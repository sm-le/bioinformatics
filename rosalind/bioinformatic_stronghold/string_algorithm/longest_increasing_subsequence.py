"""
# String algorithm 12

Given a positive integer following by a permutation of the integer,
return both longest increasing and decreasing subsequence

Example,
IN: 
    5
    5 1 4 2 3
OUT: 
    1 2 3
    5 4 2
    """

def LIS(sequence_list: list) -> str:
    """
    Find longest increasing subsequence (LIS) given a 
    list of elements

    Args:
        sequence_list: a list of elements

    Returns:
        str(increasing order)
    """

    longest_subsequence = [ [] for _ in range(len(sequence_list)) ]
    longest_subsequence[0].append(sequence_list[0])

    for idx_a in range(1, len(sequence_list)):
        for idx_b in range(idx_a):
            if int(sequence_list[idx_b]) < int(sequence_list[idx_a]) \
                and len(longest_subsequence[idx_a]) < len(longest_subsequence[idx_b]):

                longest_subsequence[idx_a] = list()
                longest_subsequence[idx_a].extend(longest_subsequence[idx_b])

        longest_subsequence[idx_a].append(sequence_list[idx_a])

    return " ".join(max(longest_subsequence, key=lambda x: len(x)))

def LDS(sequence_list: list) -> str:
    """
    Find longest decreasing subsequence (LDS) given a 
    list of elements

    Args:
        sequence_list: a list of elements

    Returns:
        str(decreasing order)
    """

    longest_subsequence = [ [] for _ in range(len(sequence_list)) ]
    longest_subsequence[0].append(sequence_list[0])

    for idx_a in range(1, len(sequence_list)):
        for idx_b in range(idx_a):
            if int(sequence_list[idx_b]) > int(sequence_list[idx_a]) \
                and len(longest_subsequence[idx_a]) < len(longest_subsequence[idx_b]):

                longest_subsequence[idx_a] = list()
                longest_subsequence[idx_a].extend(longest_subsequence[idx_b])

        longest_subsequence[idx_a].append(sequence_list[idx_a])

    return " ".join(max(longest_subsequence, key=lambda x: len(x)))

def LS(sequence_list: list) -> dict:
    """
    Find longest increasing (LIS) and decreasing subsequence (LDS) given a 
    list of elements

    Args:
        sequence_list: a list of elements

    Returns:
        dict(increasing order:str(), 
            decreasing order: str())
    """

    i_sub = LIS(sequence_list)
    d_sub = LDS(sequence_list)

    return {
        "i":i_sub,
        "d":d_sub
    }