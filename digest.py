#!/usr/bin/env python3
import re

### RESTRICTION DIGEST FUNCTION

# Inputs
seq = 'GAATTCCAAGTTCTTGTGCGCTTAAGACACAAATCCAATAAAAACTATTGTGCACACAGAATTCGCGACTTCGCGGTCTCGCTTGTTCTT' # fasta_dict[seq]
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

    ### Find all matches in forward strand
    marked_seq_fw = ""
    for match in re.finditer(rf'({pattern_left})({pattern_right})', seq):
        match_seq_fw = match.group(0)
        cut_left_fw = match.group(1) # Sequence before cut site
        cut_right_fw = match.group(2) # Sequence after cut site
        match_seq_fw_marked = (cut_left_fw + "^" + cut_right_fw) # To insert cut site symbol
    
    ### Consider reverse strand
    comp_seq = seq.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
    rv_comp_seq = comp_seq.upper()[::-1]
    print(rv_comp_seq)

    ### Find all matches in reverse strand
    marked_seq_rv = ""
    for match in re.finditer(rf'({pattern_left})({pattern_right})', rv_comp_seq):
        match_seq_rv = match.group(0)
        cut_left_rv = match.group(1) # Sequence before cut site
        cut_right_rv = match.group(2) # Sequence after cut site
        match_seq_rv_marked = (cut_left_rv + "^" + cut_right_rv) # To insert cut site symbol

    ### Replace the uncut motif with the marked version
    marked_seq_fw = seq.replace(match_seq_fw, match_seq_fw_marked)
    marked_seq_rv = rv_comp_seq.replace(match_seq_rv, match_seq_rv_marked)
    #print(marked_seq_fw)
    #print(marked_seq_rv)

    ### Split the sequence into fragments at the cut site
    fragments_fw = marked_seq_fw.split("^")
    fragments_rv = marked_seq_rv.split("^")
    #print(f'In the forward stand the fragments are {fragments_fw}')
    #print(f'In the reverse strand the fragments are {fragments_rv}')
    fragments = fragments_fw + fragments_rv
    return fragments

digest = digest(seq)
print(digest)
