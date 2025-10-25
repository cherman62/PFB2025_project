#!/usr/bin/env python3
import re

#######################################################################################################
### RESTRICTION DIGEST FUNCTION #######################################################################
#######################################################################################################

# Inputs
fasta_dict = {
    'seq1':'GAATTCCAAGTTCTTGTGCGCTTAAGACACAAATCCAATAAAAACTATTGTGCACACAGAATTCGCGACTTCGCGGTCTCGCTTGTTCTT',
    'seq2':'GAATTCCAAGTTCGATCC',
    'seq3':'AAAAAAA',
    'seq4':'ATGCGCCCGCATTT'
 } # dictionary of seq_id:sequence
motifs = ['R^AATTY', 'C^GATCC', 'GC^GGCCGC'] # list of restriction motifs, in this case for ApoI, BamHI and NotI

# Dictionary of abbreviations for all bases
bases_dict = {
    'A' : 'A', 'C' : 'C', 'G' : 'G', 'T' : 'T',
    'W' : '[AT]', 'S' : '[CG]', 'M' : '[AC]', 'K' : '[GT]',
    'R' : '[AG]', 'Y' : '[CT]',
    'B' : '[CGT]', 'D' : '[AGT]', 'H' : '[ACT]', 'V' : '[ACG]',
    'N' : '[ACGT]',
    '^' : '^' # Cut site added here to be able to convert first and split afterwards
}

### DEFINITION OF FUNCTIONS TO BE USED ################################################################
### Function to do digest 
def digest(fasta_dict, motifs):
    
    fragments_dict = {} # Empty dictionary for final output

    ### Loop over each seq in dictinory
    for seq_id, seq in fasta_dict.items():
        seq = seq.upper()
        
        ### Convert motif to a regular expression pattern
        for motif in motifs:
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
            if re.search(rf'({pattern_left})({pattern_right})', seq): # Check if ONE match exists
                #print(f'Match found for {motif} on the forward strand of {seq_id}')
                marked_seq = ""
                for match in re.finditer(rf'({pattern_left})({pattern_right})', seq): # Check ALL matches
                    match_seq  = match.group(0)
                    cut_left = match.group(1) # Sequence before cut site
                    cut_right = match.group(2) # Sequence after cut site
                    match_seq_marked = (cut_left + "^" + cut_right) # To insert cut site symbol
            else:
                #print(f'No matches found for {motif} on the forward strand of {seq_id}')
                continue      
            
            ### Replace seq with marked parts before repeating the loop
            seq = seq.replace(match_seq, match_seq_marked)

        ### Split the sequence into fragments at the cut site
        fragments = seq.split("^")
        fragments_dict[seq_id] = fragments 

    return fragments_dict

### Function to reverse complement
def rev_comp(fasta_dict):
    rev_fasta_dict = {}

    ### Loop over each seq in dictinory
    for seq_id, seq in fasta_dict.items():
        seq = seq.upper()

        ### Consider reverse strand
        seq = seq.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
        seq = seq.upper()[::-1]
        
        ### Add this to a dictionary
        rev_fasta_dict[seq_id] = seq
    return rev_fasta_dict


### BEGINNING OF ACTUAL SCRIPT USING THESE FUNCTIONS ##################################################
### Do digest
digest_dict = digest(fasta_dict, motifs) # To digest forward strand
#print(f'This is the dictionary of fragments for the forward strand {digest_dict}')
rev_fasta_dict = rev_comp(fasta_dict) # To make reverse complement
rev_digest_dict = digest(rev_fasta_dict, motifs) # To digest reverse strand
#print(f'This is the dictionary of fragments for the reverse strand {rev_digest_dict}')

### Combine dictionaries together
for key in digest_dict:
    if key in rev_digest_dict:
        digest_dict[key] += rev_digest_dict[key]
print(digest_dict)        

