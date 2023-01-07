# @sungml92
# Week 4 - Communicating with outside



######################################################
# Lecture 7.1: Communicating with the Outside Part 1 #
######################################################


'''
Reading and Writing files

1. read from file
- open a file
- open a file, but add exception when the file does not exist

2. write into a file
- 'w' to write a new file
- 'a' to append into an existing file
'''
# Reading from file
f=open('myfile','r')

try: # if you are not sure if file exists
    f=open('myfile')
except:
    print('file does not exist')

# Writing into a file
f=open('myfile','w') # to write a new file

f=open('myfile','a') # or append to existing file

'''
read content from file

1. loop over a object

2. .read() method to read a file
'''
# loop over lines in file
for line in f:
    print(line)

# basic file read
f.read() # but this does not result anything, just go to last position of a file

# change position and read again
f.seek(offset, from_what)
f.read() will begin to show line

# to read one line at a time
f.seek(0)
f.readline()

# Writing content into a file
f.write('what-to-write')

# after any file operation
f.close() # to close the file

######################################################
# Lecture 7.2: Communicating with the Outside Part 2 #
######################################################

'''
Here we will build a dictionary containing all the sequence from a fasta file

FASTA format
> id1
seq1
> id2
seq2

First we will open a file
if line is a header -> get sequence name
create new entry in dictionary
if line is not a header -> update sequence in a dictionary
if more line, continue
close file
'''

try:
    f=open("myfile.fa")
except IOError:
    print("file does not exist")

seqs = {}

for line in f:
    line = line.rstrip()

    if line[0] == ">": or # line.startswith('>')
        words=line.split()
        name = word[0][1:]
        seqs[name] = ''
    else:
        seqs[name] = seqs[name] + line
f.close()

######################################################
# Lecture 7.3: Communicating with the Outside Part 3 #
######################################################

'''
How to retrieve from dictionary

We can retrieve the key and value from dictionary using items() method

e.g)
for name,seq in seqs.items():
    print(name,seq)

Command Line Arguments
- Scripts often need to process command line arguments.
    Suppose a script that parses a FASTA file is called processfasta.py and you want to run it on a file
    > python processfasta.py myfile.fa

- The arguments of the above command are stored in the sys module's argv attribute as a list
    import sys
    print(sys.argv)

Inside the .py file
1. python path
2. comment of .py program
3. import sys
   filename=sys.argv[1] # gives filename

Parsing command line arguments with getopt
- Python's getopt module can help with processing the arguments of sys.argv
- Suppose the processfasta.py script reads a FASTA file but only stores in the dictionary the sequecnes bigger than a give length
    > processfasta.py -l 250 myfile.fa

To do this we have to tell program how to do with the argument
def usage():
    print("""

    processfasta.py [-h] [-l <length>] <filename>

    """)
Let's get into getopt

import sys
import getopt
def usage()....

o, a = getopt.getopt(sys.argv[1:], 'l:h') # o: list of optional arguments, a = list of required arguments # l:h tells l parameter expects h value
opts = {}

for k,v in o:
    opts[k] = v
if '-h' in opts.keys():
    usage(); sys.exit()

if len(a) < l:
    usage(); sys.exit("input fasta file is missing")

if '-l' in opt.keys():
    if int(opts['l']) < 0:
        print("Length of sequences should be positive"); sys.exit(0)
    seqlen=opts['-l']

Using System Environment
- When we run a script/program in the Unix environment there are standard streams recognized by a computer program
1. stdin or standard input, is stream data going into a program. Unless redirected, standard input is expected from the keyboard which started the program
2. stdout or standard output, is stream where a program writes its output data
3. stderr or standard error is another output stream used by programs to output error messages

So we can do, my_program | my_script.sh 1>program_output.txt 2>error_messages.txt

The sys module in Python provides file handles for the standard input, output and error
sys.stdin.read()
a line
another line

sys.stdout.write("Some output")

sys.stderr.write("Warning")


Interfacing with external programs
- You can call/execute an external program from within your script
- help you to automate certain task

- use call() function in the subprocess module to run an external program
    import subprocess
    subprocess.call(["ls","-l"])

- more realistically
    subprocess.call(["tophat","genome_mouse_idx","PE_reads_1.fq.gz","PE_reads_2.fq.gz"])
'''

########################
# Lecture 8: BioPython #
########################

'''
Biopython, Python based software for bioinformatics

To check installation

    import Bio
    print(Bio.__version__)

Usage Example 1: what species from unknown DNA sequence came from

my_seq.fa
-----------------
>sequence_unknown
ATCCGGT (or some sequence)
-----------------

To find out where the sequence came from, we will use "blast"

    from Bio.Blast import NCBIWWW
    fasta_string = open("myseq.fa").read()
    result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)

    to learn more, help(NCBIWWW.qblast())

To interpret result

    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)

Parsing BLAST output

    len(blast_record.alignments)

    E_VALUE_THRESH = 0.01
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print("****Alignment****")
                print('sequence',alignment.title)
                print('length:', alignment.length)
                print('e value', hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)

'''

from Bio.Blast import NCBIWWW
fasta_string = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"
result_handle = NCBIWWW.qblast("blastn","nt",fasta_string)

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print('sequence',alignment.title)
            print('length:', alignment.length)
            print('e value', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

from Bio import Seq
fasta_string2 = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"
Seq.translate(fasta_string2)
