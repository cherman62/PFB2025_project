#!/usr/bin/env python3

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

if len(sys.argv) < 2:
    print("Usage:python script.py <probe sequence>")
    sys.exit(1)

probe = sys.argv[1].upper()  # getting the probe sequence from the command line
probe_color ='yellow'

##########################################################################




############################### HTML OUTPUT FOR PROBE MATCH ###############################

def highlight_probe_html(seq, probe):
    """
    Returns the sequence as HTML with the probe highlighted in yellow and lowercase.
    """
    seq_upper = seq.upper()
    probe_upper = probe.upper()
    
    if probe_upper in seq_upper:
        start_index = seq_upper.find(probe_upper)
        end_index = start_index + len(probe_upper)
        # Keep original sequence for display, but make probe lowercase
        highlighted_seq = (
            seq[:start_index] +
            f"<span style='background-color: yellow; font-weight: bold;'>{seq[start_index:end_index].lower()}</span>" +
            seq[end_index:]
        )
        return highlighted_seq
    return None  # if no match

# Generate HTML file
with open("probe_matches.html", "w") as html_file:
    html_file.write("<html><head><title>Probe Matches</title></head><body style='background-color:white; font-family:monospace;'>\n")
    html_file.write(f"<h2>Probe Matches for '{probe}'</h2>\n")

    matches_found = False
    for sample_name, sequences in fragment_sizes.items():
        for i, seq in enumerate(sequences, start=1):
            highlighted = highlight_probe_html(seq, probe)
            if highlighted:  # only show sequences containing the probe
                html_file.write(f"<p><b>{sample_name} - Fragment {i} ({len(seq)} bp):</b> {highlighted}</p>\n")
                matches_found = True

    if not matches_found:
        html_file.write("<p>No sequences matched the probe.</p>\n")

    html_file.write("</body></html>")

print("HTML file 'probe_matches.html' created with probe-highlighted sequences.")

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
sample_color = "cyan"
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
        #if the probe is found in the fragment, color that band yellow
        if probe in seq: 
            band_color = probe_color
            print(f"Probe '{probe}'detected in {sample_name} fragment ({len(seq)} bp)") ##printing to screen the probe "seq", and where it is detected 
        else: 
            band_color = sample_color

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
        ax.hlines(y, fragment_x-0.2, fragment_x+0.2, color=band_color, linewidth=6)
        ax.text(fragment_x+0.3, y, f"{len(seq)} bp", color='white', va='center', fontsize=8)

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
with open("dna_fragments_lengths.txt", "w") as f:
    for sample_name, sequences in fragment_sizes.items():
        f.write(f"Sample: {sample_name}\n")
        for i, seq in enumerate(sequences, start=1):
            length = len(seq)
            f.write(f"  Fragment {i}: Length = {length} bp\n")
            f.write(f"    Sequence: {seq}\n")
        f.write("\n")  # extra line between samples
#############################################################################



