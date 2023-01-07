"""
# Dynamic programming 1

Given positive integers n≤40 and k≤5, return the total number 
of rabbit pairs that will be present after n months, k

Example,
IN: 5 3
OUT: 19
"""

def fibonacci_v1(month:int, offspring:int) -> int:
    """Fibonacci sequence with recurrence relation

    Args:
        month: month
        offspring: how many offspring pairs does a pair of rabbit produce?
    
    Returns:
        The total number of rabbit pairs after n months
    """

    if month <= 2:
        return 1
    
    return fibonacci_v1(month-1,offspring) + (fibonacci_v1(month-2,offspring) * offspring)

def fibonacci_v2(month:int, offspring:int) -> int:
    """Fibonacci sequence with recurrence relation

    Args:
        month: month
        offspring: how many offspring pairs does a pair of rabbit produce?
    
    Returns:
        The total number of rabbit pairs after n months
    """

    
    return 1 if month <= 2 else \
        fibonacci_v2(month-1,offspring) + (fibonacci_v2(month-2,offspring) * offspring)
    