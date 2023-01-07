"""
# introduction 1

Given two positive integers a, b <= 1000, return hypotenuse of the right triangle

Example,
IN: 3 5
OUT: 34
"""

def hypotenuse_square(one_side:int, other_side:int) -> int:
    """Calculate square of hypotenuse of the right triangle

    Args:
        one_side: length of one side of the triangle
        other_side: length of the other side of the triangle

    Returns:
        length of hypotenuse
    """

    hypotenuse = (one_side**2) + (other_side**2)

    return hypotenuse