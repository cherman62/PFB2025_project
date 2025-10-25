#!/usr/bin/env python3
import re

### RESTRICTION DIGEST FUNCTION

# Inputs
seq = 'GAATTCCAAGTTCTTGTGCGCACACAAATCCAATAAAAACTATTGTGCACACAGAATTCGCGACTTCGCGGTCTCGCTTGTTCTT' # fasta_dict[seq]
motif = 'R^AATTY' # enz_dict[enzyme]

# Dictionary of abbreviations for all bases
bases_dict = {
    'A' : 'A', 'C' : 'C', 'G' : 'G', 'T' : 'T',
    'W' : '[AT]', 'S' : '[CG]', 'M' : '[AC]', 'K' : '[GT]',
    'R' : '[AG]', 'Y' : '[CT]',
    'B' : '[CGT]', 'D' : '[AGT]', 'H' : '[ACT]', 'V' : '[ACG]',
    'N' : '[ACGT]',
    '^' : '^' # Cut site added here to be able to convert first and split afterwards
}
       
# Function to do restriction digest
def digest(seq):
    # for line in seq:
    #     line = line.rstrip () # to remove the empty character

    seq = seq.upper()

    ### Convert motif to a regular expression pattern
    pattern = ''
    for base in motif:
        if base in bases_dict:
            pattern += bases_dict[base] # Replace degenerate bases with regex equivalent
        if base not in bases_dict:
            print(f'Error: Check spelling of motif')
            exit(1)       
    
    ### Split regex pattern into two at cut site symbol for matching
    pattern_left, pattern_right = pattern.split("^")        

    ### Find all matches
    marked_seq = ""
    for match in re.finditer(rf'({pattern_left})({pattern_right})', seq):
        match_seq = match.group(0)
        cut_left = match.group(1) # Sequence before cut site
        cut_right = match.group(2) # Sequence after cut site
        match_seq_marked = (cut_left + "^" + cut_right) # To insert cut site symbol
    
    ### Replace the uncut motif with the marked version
    marked_seq = seq.replace(match_seq, match_seq_marked)
    #print(marked_seq)

    ### Split the sequence into fragments at the cut site
    fragments = marked_seq.split("^")
    return fragments

digest = digest(seq)
print(digest)
