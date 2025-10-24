#!/usr/bin/env/ python3 

import Bio
import Bio.Seq
from Bio import SeqIO

def read_fasta(filename):

    seqID_seq = {}

    for seq_record in SeqIO.parse(filename, 'fasta'):
        seqID_seq[seq_record.id] = seq_record.seq 

    return seqID_seq

print(read_fasta('multiple_seq_test.fasta'))