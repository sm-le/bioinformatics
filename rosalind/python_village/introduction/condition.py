"""
# introduction 3

Given two integers a and b, sum all odd number between two integer

Example,
IN: 100 200
OUT: 7500
"""

def oddnumber_sum(a:int, b:int) -> int:
    """Calculate square of hypotenuse of the right triangle

    Args:
        a: start number
        b: end number

    Returns:
        sum of odd numbers between range of a,b
    """

    odd_sum = sum([i for i in range(a, b+1) if i % 2 == 1])

    return odd_sum