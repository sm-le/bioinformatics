"""
# Dynamic programming 2

Given positive integers n≤100 and m≤20, return the total number 
of rabbit pairs that will remain after n months with life expectancy 
of m

Example,
IN: 6 3
OUT: 4
"""

def mortal_fibonacci_v1(month:int, life_span:int) -> int:
    """Fibonacci sequence with recurrence relation and mortality

    Args:
        month: month
        life_span: how long does a pair of rabbit survives
    
    Returns:
        The remaining number of rabbit pairs after n months
    """
    if month <= 2:
        return 1 if month > 0 else 0

    elif month <= life_span:
        return mortal_fibonacci_v1(month-1,life_span) + mortal_fibonacci_v1(month-2,life_span)
    
    elif month == life_span + 1:
        return mortal_fibonacci_v1(month-1,life_span) + mortal_fibonacci_v1(month-2,life_span) - 1

    return mortal_fibonacci_v1(month-1,life_span) + mortal_fibonacci_v1(month-2,life_span) - mortal_fibonacci_v1(month-(life_span + 1), life_span)

def mortal_fibonacci_v2(month:int, life_span:int) -> int:
    """Fibonacci sequence with recurrence relation and mortality

    Args:
        month: month
        life_span: how long does a pair of rabbit survives
    
    Returns:
        The total number of rabbit pairs after n months
    """
    my_bunnies = list()
    current_month = 1
    try:
        while current_month <= month:
            if current_month <= 2:
                my_bunnies.append(1)
            elif current_month <= life_span:
                my_bunnies.append(my_bunnies[-1]+my_bunnies[-2])
            elif current_month == life_span + 1:
                my_bunnies.append(my_bunnies[-1]+my_bunnies[-2]-1)
            else:
                my_bunnies.append(my_bunnies[-1]+my_bunnies[-2]-my_bunnies[-(life_span+1)])
            
            current_month += 1

        return my_bunnies[-1]
        
    except Exception as e:
        print(e)