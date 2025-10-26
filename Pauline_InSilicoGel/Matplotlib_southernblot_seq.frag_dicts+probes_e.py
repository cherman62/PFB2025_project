#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math

#DATA
fragment_sizes = {"Seq1": ["ATGCGT" * 500 , "GATTACA" * 300 , "TTAGGC" * 270 , "CCGTA" *100, "AT" *72],
                  "Seq2": ["GCGT" * 300 , "TTGAC" * 312 , "ATGC" * 120 , "CG" * 75]}

ladder = {'100bp': 100, '200bp': 200, '500bp': 500,
          '1000bp': 1000, '2000bp': 2000, '3000bp': 3000}

probes = {"Probe1" : "ATGCGT"*100,
          "Probe2" : "GATTACA"*50 , 
          "Probe3" : "TTGAC"*25 , 
          "Probe4" : "ATGC"*200 , 
          "Probe5" : "GCGT" *150}

probe_colors = ["yellow" , "lime" , "magenta" , "cyan" , "orange", "teal" , "blue", "coral" , "darkgreen" , "tan", "sienna", "plum", "lavendar"]



############################################################# FUNCTIONS #####################################################################

#Simulate migration distance (smaller fragments migrate farther)- plotted on a logarithmic scale, as that is how typical DNA agarose gels run 
def migration_distance(size):
    return 120 - (math.log10(size) * 25)  # inverted log spacing

#Highlighting the probe in HTML format 
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

#Setting up and generating HTML file 
def generate_html(filename, fragment_sizes, probes, probe_colors): 
    matches_found = False
    with open("probe_matches_e.html", "w") as html_file:
        html_file.write("<html><head><title>Probe Matches</title></head><body style='background-color:white; font-family:monospace;'>\n")
        html_file.write(f"<h2>Probe Matches with Tm Values</h2>\n")

        # ----------- Add probe legend with Tm values -----------
        html_file.write("<h3>Probe Legend</h3>\n<ul>\n")
        for i, (probe_name, probe_seq) in enumerate(probes.items()):
            color = probe_colors[i % len(probe_colors)]
            tm = calculate_tm(probe_seq)
            html_file.write(
                f"<li><b>{probe_name}</b>: "
                f"<span style='background-color:{color};'>{probe_seq}</span> "
                f"({len(probe_seq)} bp, Tm = {tm}C)</li>\n")
        html_file.write("</ul><hr>\n")

        # -------Highlight probe matches in sequences ------------
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

#DNA ladder in gel
def plot_ladder(ax, ladder, x_pos, color ="red"):
    for name, size in ladder.items():
        y = migration_distance(size)
        ax.hlines(y, x_pos-0.2, x_pos+0.2, color=color, linewidth=3)
        ax.text(x_pos-0.5, y, name, color='white', va='center')
    ax.text(x_pos, +20, "DNA Ladder", color=color, ha='center', fontsize=12, fontweight='bold')

#Sample lane in gel 
def plot_sample_lane(ax, fragment_sizes, probes, probe_colors, x_pos, lane_lable="Test Sample"):
    sample_color = "white"
    for sample_name, seq_list in fragment_sizes.items(): 
        for seq in seq_list:
            y = migration_distance(len(seq))
            band_color = sample_color
            for i, probe_seq in enumerate(probes.values()):
                if probe_seq in seq:
                    band_color = probe_colors[i % len(probe_colors)]
                    print(f"Probe '{probe_seq}' detected in {sample_name} fragment ({len(seq)} bp)")
                    break  # stop at first matching probe
            ax.hlines(y, x_pos-0.2, x_pos+0.2, color=band_color, linewidth=2)
        ax.text(x_pos, +20, "Test Sample", color=sample_color, ha='center', fontsize=12, fontweight='bold')

#Saving fragment lengths to a file(.txt)
def save_fragment_lengths(filename, fragment_sizes): 
    with open("dna_fragments_lengths_d.txt", "w") as f:
        for sample_name, sequences in fragment_sizes.items():
            f.write(f"Sample: {sample_name}\n")
            for i, seq in enumerate(sequences, start=1):
                f.write(f"  Fragment {i}: Length = {len(seq)} bp\n Sequence: {seq}\n")
            f.write("\n")  # extra line between samples

#Melting temperature calculation (Wallace Rule)
def calculate_tm(seq): 
    seq = seq.upper()
    a = seq.count ("A")
    t = seq.count ("T")
    g = seq.count ("G")
    c = seq.count ("C")
    total = a + t + g + c 
    if total < 14: 
        tm = 2*(a+t) + 4*(g+c) # using Wallace rule for shorter (less than 14 bp) oligos 
    else:
        tm = 64.9 + 41*(g+c -16.4)/total
    return round (tm, 1)


############################################## MAIN SCRIPT ####################################################

#Generate HTML with Tm values
generate_html("probe_matches_e.html", fragment_sizes, probes, probe_colors)

# Gel plot setup
fig, ax = plt.subplots(figsize=(6, 10))
ax.set_facecolor('black')
ax.set_ylim(0,130)
ax.set_xlim(0,2)#adjusted to fir both lanes 
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])

#Plot ladder 
ladder_x = 0.5
plot_ladder(ax, ladder, ladder_x)

#Plot sample lane 
lane_x = 1.5
plot_sample_lane(ax, fragment_sizes, probes, probe_colors, lane_x)

#Title of plot 
title_x = (ladder_x + lane_x) / 2
ax.text(title_x, 5, "In Silico DNA Gel", color='white', ha='center', fontsize=16, fontweight='bold')

#Saving plot
plt.savefig("in_silico_test_gel_fixed_e.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

#Save fragment lengths
save_fragment_lengths("dna_fragments_lengths_d.txt", fragment_sizes)







# if sequences are coming from FASTA file format, sequences with multiple fragments 

#Replace harcoded fragment sizes
# # Load fragment sizes from FASTA
# fragment_sizes = fasta_to_fragment_sizes_grouped("all_fragments.fasta")


#Add FASTA parser at the top
# Function to read fragments from a FASTA file and group by sample
# def fasta_to_fragment_sizes_grouped(fasta_file):
#     fragment_sizes = {}
#     with open(fasta_file, "r") as f:
#         current_header = None
#         seq_lines = []
#         for line in f:
#             line = line.strip()
#             if line.startswith(">"):
#                 if current_header:
#                     # Extract sample name (before first underscore)
#                     sample_name = current_header.split("_")[0]
#                     fragment_sizes.setdefault(sample_name, []).append("".join(seq_lines))
#                 current_header = line[1:]  # remove ">"
#                 seq_lines = []
#             else:
#                 seq_lines.append(line)
#         # Add last sequence
#         if current_header:
#             sample_name = current_header.split("_")[0]
#             fragment_sizes.setdefault(sample_name, []).append("".join(seq_lines))
#     return fragment_sizes





# if probes come in as a FASTA format

# replace harcoded probes: 
# # Load probes from FASTA
# probes = fasta_to_probes("probes.fasta")

# make a function that takes probes from a FASTA file 
# def fasta_to_probes(fasta_file):
#     """
#     Reads a FASTA file containing probes and returns a dictionary:
#     {header: sequence}
#     """
#     probes = {}
#     with open(fasta_file, "r") as f:
#         current_header = None
#         seq_lines = []
#         for line in f:
#             line = line.strip()
#             if line.startswith(">"):
#                 if current_header:
#                     probes[current_header] = "".join(seq_lines)
#                 current_header = line[1:]  # remove ">"
#                 seq_lines = []
#             else:
#                 seq_lines.append(line)
#         # Add last probe
#         if current_header:
#             probes[current_header] = "".join(seq_lines)
#     return probes
