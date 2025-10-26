#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math
import sys

############### DICT OF FRAGMENT AND THEIR SEQUENCE LENGTH #################

fragment_sizes = {"Seq1": ["ATGCGT" * 500 , "GATTACA" * 300 , "TTAGGC" * 270 , "CCGTA" *100, "AT" *72],
                  "Seq2": ["GCGT" * 300 , "TTGAC" * 312 , "ATGC" * 120 , "CG" * 75]}
ladder = {'100bp': 100, '200bp': 200, '500bp': 500,
          '1000bp': 1000, '2000bp': 2000, '3000bp': 3000}

def migration_distance(size):
   #"""Simulate migration distance (smaller fragments migrate farther)."""
    return 120 - (math.log10(size) * 25)  # inverted log spacing
############################################################################





########################### ADD PROBE HANDLING ###########################

probes = {"Probe1" : "ATGCGT",
          "Probe2" : "GATTACA" , 
          "Probe3" : "TTGAC"}

#define colors
probe_colors = ["yellow" , "lime" , "magenta" , "cyan" , "orange"]


##########################################################################




############################### HTML OUTPUT FOR PROBE MATCH ###############################

def highlight_probe_html(seq, probe_dict, colors):
    """
    Returns the sequence as HTML with all probes highlighted in different colors.
    """
    seq_upper = seq.upper()
    highlighted_seq = seq  # start with original sequence

    for i, (probe_name, probe_seq) in enumerate(probe_dict.items()):  # <- use probe_dict
        probe_upper = probe_seq.upper()
        if probe_upper in seq_upper:
            start_index = seq_upper.find(probe_upper)
            end_index = start_index + len(probe_upper)
            color = colors[i]
            highlighted_seq = (
                highlighted_seq[:start_index] +
                f"<span style='background-color: {color}; font-weight: bold;'>{highlighted_seq[start_index:end_index].lower()}</span>" +
                highlighted_seq[end_index:])
            # update seq_upper to avoid overlapping replacements
            seq_upper = highlighted_seq.upper()
    return highlighted_seq
   

# Generate HTML file
with open("probe_matches_d.html", "w") as html_file:
    html_file.write("<html><head><title>Probe Matches</title></head><body style='background-color:white; font-family:monospace;'>\n")
    html_file.write(f"<h2>Probe Matches</h2>\n")

    matches_found = False
    for sample_name, sequences in fragment_sizes.items():
        for i, seq in enumerate(sequences, start=1):
            highlighted = highlight_probe_html(seq, probes, probe_colors)
            if highlighted != seq:  # if any probe was highlighted
                html_file.write(f"<p><b>{sample_name} - Fragment {i} ({len(seq)} bp):</b> {highlighted}</p>\n")
                matches_found = True

    if not matches_found:
        html_file.write("<p>No sequences matched any probe.</p>\n")

    html_file.write("</body></html>")

print("HTML file and legend for probes in gel 'probe_matches_d.html' created with probe-highlighted sequences.")

############################################################################





############################### PLOT SETUP ################################

# Plot setup
fig, ax = plt.subplots(figsize=(6, 10))
ax.set_facecolor('black')
ax.set_ylim(0,130)
ax.set_xlim(0,2)#adjusted to fir both lanes 
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])

ladder_color = "red"
sample_color = "white"
############################################################################





############################# PLOT LADDER ##################################
ladder_x = 0.5
for name, size in ladder.items():
    y = migration_distance(size)
    ax.hlines(y, ladder_x-0.2, ladder_x+0.2, color=ladder_color, linewidth=3)
    ax.text(ladder_x-0.5, y, name, color='white', va='center')

#DNA ladder label below the ladder lane 
ax.text(ladder_x, +20, "DNA Ladder", color=ladder_color, ha='center', fontsize=12, fontweight='bold')
############################################################################




########################## PLOT SAMPLE LANE #################################
lane_gap = 1.0  # Gap between lanes for different samples
lane_index = 1  # Start lane index (after ladder lanes)

for sample_name, seq_list in fragment_sizes.items():
    fragment_x = ladder_x + lane_gap * lane_index  # Calculate the X position for the lane

    #Now I want to sort fragment by size to make font calculation easier (thus if there are bands clos in size, I want the font to be smaller)
    prev_label_y=None
    min_gap=4 #minimal vertical gap between labels 

    for seq in fragment_sizes[sample_name]:
        y = migration_distance(len(seq))


        #######>>>> Probe Match Coloring <<<<<######
        # Determine band color based on probes
        band_color = sample_color  # default color
        for i, probe_seq in enumerate(probes.values()):
            if probe_seq in seq:
                band_color = probe_colors[i % len(probe_colors)]
                print(f"Probe '{probe_seq}' detected in {sample_name} fragment ({len(seq)} bp)")
                break  # stop at first matching probe

    # Plot all fragments for this sample in the corresponding lane
        ax.hlines(y, fragment_x-0.2, fragment_x+0.2, color=band_color, linewidth=2)

        prev_y = y

    # Label the lane with the sample name
    ax.text(fragment_x, +20, "Test Sample", color=sample_color, ha='center', fontsize=12, fontweight='bold')
#############################################################################





################################# TITLE ######################################
title_x = (ladder_x + fragment_x) /2 #this will calculate the center point between the two lane to put the label 
ax.text(title_x, +5, "In Silico DNA Gel", color='white', ha = 'center', fontsize=16, fontweight='bold')
#############################################################################



############################# PLOT GENERATION ###############################
plt.savefig("in_silico_test_gel_fixed.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()
#############################################################################



######################### PRINTING FRAG SEQS AND LENGTHS ###################
with open("dna_fragments_lengths_d.txt", "w") as f:
    for sample_name, sequences in fragment_sizes.items():
        f.write(f"Sample: {sample_name}\n")
        for i, seq in enumerate(sequences, start=1):
            length = len(seq)
            f.write(f"  Fragment {i}: Length = {length} bp\n")
            f.write(f"    Sequence: {seq}\n")
        f.write("\n")  # extra line between samples
#############################################################################



