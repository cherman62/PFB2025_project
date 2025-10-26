#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math

####### HERE I WANT EACH SEQUENCE LENGTH TO DETERMINE THE MIGRATION DISTANCE 
fragment_sizes = {"Seq1": ["ATGCGT" * 500 , "GATTACA" * 300 , "TTAGGC" * 270 , "CCGTA" *100, "AT" *72],
                  "Seq2": ["GCGT" * 300 , "TTGAC" * 312 , "ATGC" * 120 , "CG" * 75]}
ladder = {'100bp': 100, '200bp': 200, '500bp': 500,
          '1000bp': 1000, '2000bp': 2000, '3000bp': 3000}

def migration_distance(size):
   #"""Simulate migration distance (smaller fragments migrate farther)."""
    return 120 - (math.log10(size) * 25)  # inverted log spacing

###########################################################################

# Plot setup
fig, ax = plt.subplots(figsize=(6, 10))
ax.set_facecolor('black')
ax.set_ylim(0,130)
ax.set_xlim(0,2)#adjusted to fir both lanes 
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])

ladder_color = "red"
sample_color = "cyan"

############################################################################

# Plot Ladder
ladder_x = 0.5
for name, size in ladder.items():
    y = migration_distance(size)
    ax.hlines(y, ladder_x-0.2, ladder_x+0.2, color=ladder_color, linewidth=3)
    ax.text(ladder_x-0.5, y, name, color='white', va='center')

#DNA ladder label below the ladder lane 
ax.text(ladder_x, +20, "DNA Ladder", color=ladder_color, ha='center', fontsize=12, fontweight='bold')

############################################################################

# Plot Sample Lane 
lane_gap = 1.0  # Gap between lanes for different samples
lane_index = 1  # Start lane index (after ladder lanes

for sample_name, seq_list in fragment_sizes.items():
    fragment_x = ladder_x + lane_gap * lane_index  # Calculate the X position for the lane

    #Now I want to sort fragment by size to make font calculation easier (thus if there are bands clos in size, I want the font to be smaller)
    prev_label_y=None
    min_gap=4 #minimal vertical gap between labels 

    for seq in fragment_sizes[sample_name]:
        y = migration_distance(len(seq))
        #adjust font size based on proximity to previous band 
        if prev_label_y is None:
            label_y=y
        else:
            distance = prev_label_y - y 
            if distance < min_gap:
                label_y = prev_label_y - min_gap
            else:
                label_y = y
              

    # Plot all fragments for this sample in the corresponding lane
        ax.hlines(y, fragment_x-0.2, fragment_x+0.2, color=sample_color, linewidth=3)
        ax.text(fragment_x+0.3, y, f"{len(seq)} bp", color='white', va='center', fontsize=8)

        prev_y = y

# Label the lane with the sample name
    ax.text(fragment_x, +20, "Test Sample", color=sample_color, ha='center', fontsize=12, fontweight='bold')

#############################################################################


# Title
title_x = (ladder_x + fragment_x) /2 #this will calculate the center point between the two lane to put the label 
ax.text(title_x, +5, "In Silico DNA Gel", color='white', ha = 'center', fontsize=16, fontweight='bold')

#############################################################################

plt.savefig("in_silico_test_gel_fixed.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

#############################################################################
#Becuase some bands overlap, labels are also going to overlap so to be safe I want to put all fragment sequences into a txt file and have their corresponding length. 

with open("dna_fragments_lengths.txt", "w") as f:
    for sample_name, sequences in fragment_sizes.items():
        f.write(f"Sample: {sample_name}\n")
        for i, seq in enumerate(sequences, start=1):
            length = len(seq)
            f.write(f"  Fragment {i}: Length = {length} bp\n")
            f.write(f"    Sequence: {seq}\n")
        f.write("\n")  # extra line between samples