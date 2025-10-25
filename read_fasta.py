#!/usr/bin/env/ python3 

def read_fasta(fasta_sequence):

    fasta_file_split = fasta_sequence.split('\n')
    filtered_list = [real_string for real_string in fasta_file_split if real_string.strip()]

    seq_ID = []
    seq = []
    growing_str = ''

    for list_element in filtered_list:
        if list_element.startswith('>'):
            if growing_str:
                seq.append(growing_str)
                growing_str = ''
            seq_ID.append(list_element)
        else:
            growing_str += list_element
    seq.append(growing_str)

    fasta_dict = dict(zip(seq_ID,seq))

    return fasta_dict, len(fasta_dict)

example_fasta_file = '''
>seq1
AAGAGCAGCTCGCGCTAATGTGATAGATGGCGGTAAAGTAAATGTCCTATGGGCCACCAATTATGGTGTATGAGTGAATCTCTGGTCCGAGATTCA
CTGAGTAACTGCTGTACACAGTAGTAACACGTGGAGATCCCATAAGCTTCACGTGTGGTCCAATAAAACACTCCGTTGGTCAAC
>seq2
GCCACAGAGCCTAGGACCCCAACCTAACCTAACCTAACCTAACCTACAGTTTGATCTTAACCATGAGGCTGAGAAGCGATGTCCTGACCGGCCTGT
CCTAACCGCCCTGACCTAACCGGCTTGACCTAACCGCCCTGACCTAACCAGGCTAACCTAACCAAACCGTGAAAAAAGGAATCT
>seq3
ATGAAAGTTACATAAAGACTATTCGATGCATAAATAGTTCAGTTTTGAAAACTTACATTTTGTTAAAGTCAGGTACTTGTGTATAATATCAACTAA
AT
>seq4
ATGCTAACCAAAGTTTCAGTTCGGACGTGTCGATGAGCGACGCTCAAAAAGGAAACAACATGCCAAATAGAAACGATCAATTCGGCGATGGAAATC
AGAACAACGATCAGTTTGGAAATCAAAATAGAAATAACGGGAACGATCAGTTTAATAACATGATGCAGAATAAAGGGAATAATCAATTTAATCCAG
GTAATCAGAACAGAGGT
'''

print(read_fasta(example_fasta_file))