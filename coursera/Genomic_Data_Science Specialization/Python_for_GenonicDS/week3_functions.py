# @sungml92
# Week 3 - Functions



#######################
# Lecture 5 Functions #
#######################
'''
Lecture 5.1 Functions part 1

A functions is a part of program. It takes a list of argument values,
    performs a computation with those values, and returns a single result

There are given function like len() or print()

Why write functions?

 1. Reusability
        allow us to reuse code instead of rewriting it
 2. Abstraction
        allow us to conceive program as a sequence of sub-steps.

Some useful application of function in bioinformatics

 1. A function that computes GC contents

 2. A function that checks in-frame stop codon

 3. A function that reverse complement a DNA sequence

General structure of function follow

    def funtion_name(input arguments):
        "string documenting this function"

        function code block

        return output # return ouput we want

'''

# GC content function

def gc(dna):
    "This function computes GC percentage of a dna sequence"
    dna = dna.lower()
    nbases = dna.count('n')
    gcpercent = (dna.count('c') + dna.count('g') * 100.0) / (len(dna) - nbases)

    return gcpercent

gc('AAGCTTGGAA')
help(gc)

'''
Scope of variable declaration

 If you declare a variable inside the function it will only exist in that function
    For example, nbases in gc function cannot be accessed in main prompt

Lecture 5.2 Function part 2

Boolean function
 are functions that return a True or False
'''
# Example Boolean function, check in-frame stop codon

def dna_has_stop(dna):
    "This function checks if given dna sequence has in-frame stop codons"
    stop_codon_found = False
    stop_codons = ['tga','taa','tag']
    dna = dna.lower()
    for i in range(0,len(dna),3):
        codon = dna[i:i+3]
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found

'''
Defining function default values
 Suppose the has_stop_codon function also accepts a frame argument to look at frame of interest (0 or 1 or 2)
'''
def dna_has_stop(dna,frame=0):
    "This function checks if given dna sequence has in-frame stop codons"
    stop_codon_found = False
    stop_codons = ['tga','taa','tag']
    dna = dna.lower()
    for i in range(frame,len(dna),3):
        codon = dna[i:i+3]
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found

def dna_has_stop_auto(dna):
    "This function checks if given dna sequence has in-frame stop codons"
    "This function automatically detect frame that has stop codon"

    stop_codon_pos = []
    stop_codon_frame = []
    stop_codons = ['tga','taa','tag']
    dna = dna.lower()
    for frame in range(3):
        for i in range(frame,len(dna),3):
            codon = dna[i:i+3]
            if codon in stop_codons:
                stop_codon_pos.append(i)
                stop_codon_frame.append(frame)

    for a,b in zip(stop_codon_pos, stop_codon_frame):
        print("stop codon found in %d at %d' frame \n" % (a,b))
    return

dna = 'atgagcggccggttgatgatga'
dna_has_stop_auto(dna)

'''
Lecture 5.3: Function Part 3

Passing Function arguments

 1. pass arguments by position
        has_stop_codon(dna,1)

 2. pass arguments by name
        has_stop_codon(frame=0,dna=seq)

 3. mixed either style, but named argument must be put after positioned arguments
        has_stop_codon(seq,frame=2)

A More complex function example
 Write a function to reverse complement a DNA sequence

    def reversecomplement(dna):
        seq = reverse_string(seq)
        seq = complement(seq)

        return seq

Reverse a string

 1. Regular slice
    dna = "AGT"
    dna[0:3]
 2. Extended slices
    dna[0:3:2]
 3. dna[::-1] # backward step

List comprehensions
 - provide a concise way to create lists
 - Common applications are to make new lists where each element is the results of some operation

 So

 new_list = [operation(i) for i in old_list if filter(i)]
 equals to
 new_list = []
 for i in old_list:
     if filter(i):
         new_list.append(operation(i))
'''

# Reverse complement example
def reversecomplement(dna):
    "This function returns the reverse complement of the dna string"

    basecomplement = {'a':'t','c':'g','g':'c','t':'a','n':'n'}
    dna = dna.lower()

    revdna = dna[::-1]
    letters = list(revdna)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


    return seq

'''
Split and Join
 are method of the string object

 1. split() returns a list of all the words in the string

 2. join() returns a string in which the string elements were joined from a list
'''

##################################
# Lecture 6 Modules and Packages #
##################################

'''
What are Modules?
 - Modules in python are simply python files with the .py extension, which contain definitions of functions, or vars
 - Grouping related code into a module makes the code easier to understand and use

 We could put all previous functions in a file named dnautil and import the dnautil.py file by
    import dnautil

The sys.path variable
 you can use the sys.path variable from the sys built in module to check the list of all directories
 e.g) import sys
      sys.path

 If the sys.path doesn't contain the dir, you can extend it
    sys.path.append(path-to-dir)

Importing name from a module
 - you can import all names that a module defines with the following statement
    from dnautil import *

Packages
 packages group multiple modules under one name, by using 'dotted module name'.
  For example, the module name A,B designates a submodule named B in a package named A.

 Each package in python is a directory which must contain a special file called __init__.py
'''

import time
start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))

dir(time)
