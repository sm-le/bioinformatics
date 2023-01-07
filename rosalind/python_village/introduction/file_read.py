"""
# introduction 4

Given a file, read even numbered line

Example,
IN:
    yes
    I can do
    it
    you?
OUT: 
    I can do
    you?
"""

def file_reader_even(file_path:str) -> str:
    """Read even numbered row from file

    Args:
        file_path: path to file
    
    Returns:
        even numbered row
    """
    
    my_line = """"""

    with open(file_path, "r") as f:
        for row_num, line in enumerate(f):
            if row_num % 2 == 1:
                my_line += f"{line}"

    return my_line
