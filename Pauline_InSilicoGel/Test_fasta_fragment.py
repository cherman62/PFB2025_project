#!/usr/bin/env python3

#This is for a single dictionary of fragments!! 

# make a dictionary of synthetic dna fragments
test_fragments = {"Fragment_1_250bp" :"ATGCGTACGATGACCTGAGCGTTACCGATGGCTTACGCTAGGCTGATCGTAGCTGACGTACGTGACGATCGTAGCGTACGTAGCTAGCGTACGATGCTAGCGATGCGTAGCTGACGTAGCGTACGTAGCTAGCGTACGTAGCGTACGATGCTAGCGATGCGTACGTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTACGTAGCGA", 
                  "Fragment_2_500bp" :"GATCGTACGATGCTAGCGTACGTAGCGTACGATCGTAGCTAGCGTAGCTAGCGTACGTAGCGTACGATGCTAGCGATGCGTACGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTACGTAGCGTACGATGCTAGCGATGCGTAGCGTACGTAGCGTACGTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTACGATGCTAGCGATGCGTAGCTAGCGTACGTAGCTAGCGTAGCGTACGATGCTAGCGA",
                  "Fragment_3_800bp" :"ATCGTACGATGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTACGTAGCTAGCGTACGTAGCGTAGCTAGCGTACGATGCTAGCGATGCGTAGCTAGCGTAGCTAGCGTACGTAGCGTACGATGCTAGCGATGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTACGTAGCTAGCGTACGATGCTAGCGATGCGTAGCGTACGTAGCGTACGTAGCTAGCGTAGCTAGCGTAGCGTACGATGCTAGCGATGCGTACGTAGCGTACGTAGCGTACGTAGCGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTAGC" , 
                  "Fragment_4_1200bp" :("GCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTACG" * 12)[:1200] , 
                  "Fragment_5_2000bp" : ("ATCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGT" * 20)[:2000] , 
                  "Fragment_6_3000bp" : ("GATCGTACGATGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGCGTAGCTAGCGTAGCTAGCGTAGCTAGCGTAGC" * 30)[:3000]}

#now I wan to make a function that wraps lines at 80 characters (typically FASTA format standard)
#below is the function to wrap lines at 80 characters
def wrap_sequence(seq,width = 80):
    result = "" #starting with an empty string 
    for i in range (0, len(seq), width): # going through the sequence in 80b steps 
        result  += seq[i:i+width] + "\n" #take each 8- bp slice and add a newline 
    return result #return the sequence in new format 

#now I want to create FASTA content and here I am building the text that goes into the FASTA file. 
fasta_content = ""  
for name, seq in test_fragments.items():# loop through each fragment in the dictionary 
        fasta_content += f">{name}\n" + wrap_sequence(seq) # add the wrapped DNA sequence and the FASTA header line (ie., >Fragment_1_250bp)

#Save the FASTA file and define the name of the file 
filename="restriction_digest_fragments_test.fasta"
with open(filename, "w") as fh:
     fh.write(fasta_content)

#Calculare total sequence length (by suming the length of each sequence)
total_bp = 0
for seq in test_fragments.values(): 
     total_bp = total_bp + len(seq)

#print out summary information 
print("FASTA file successfully created")
print("File name:" , filename)
print("Total fragments:", len(test_fragments))
print("Total sequence length(bp):" , total_bp)