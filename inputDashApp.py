#!/usr/bin/env python

from dash import Dash, html, dcc, Input, Output, State, ALL, no_update
import re
from read_fasta import read_fasta
from read_re import re_match
from digest import *
import matplotlib.pyplot as plt
import math
from insilico_gel import *
import io 
import base64

import os

#Initialize the app
app = Dash()

# Counter for dynamically created enzyme input boxes
counter = {"n": 1}  # start from 1 because the first box has index 0

#App layout
app.layout = html.Div([
    html.Div([                                                          # Header
    html.Header("In silico Southern blotting tool", 
                 style={"font-size": "30px", "textAlign": "center"})],
     style={"marginBottom": "20px"}),
     html.Div([                                                          # Descripton
    html.Label("Southern blot is a technique used to detect target DNA. The in silico southern blot tool facilitates the user's southern blot experiment by visualizing probe binding to DNA sample cleaved by restriction enzyme(s) of choice.", 
                 style={"font-size": "20px"}),
    html.Label("To start, please enter a FASTA and probe sequence(s), and up to two restriction enzymes. Restriction enzyme names are case-sensitive.", 
                 style={"font-size": "20px"}),
    html.Label("The output is a virtual gel of digested DNA sample and displays color-coded bands indicating locations the probe(s) are expected to bind.", 
                 style={"font-size": "20px"})],
     style={"marginBottom": "30px"}),

    html.Div([                                                          # Fasta input title  
     html.Label("1) Paste one or more DNA sequences in FASTA format", htmlFor='fasta_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Div([                                                          # Fasta input box 
     dcc.Textarea(id='fasta_input', rows=20, cols=90, value='', required=True,
                placeholder=f'>sequence name\nACGTACGT...\n>sequence2 name\nACGTACGT...\n',)],
                style={"marginBottom": "5px"}),
#######################################
    html.Div(dcc.Upload(id= 'Upload button', contents='', children = html.Div(html.A('Select file',
                 style={"font-size": "20px", "marginBottom": "20px", 'marginLeft': '0px', "color": "#007bff", "textDecoration": "none"}))),
                 style={"marginBottom": "20px"}),


######################################

    html.Div([                                              # Probe Fasta input title  
     html.Label("2) Paste one or more Probe sequences in FASTA format", htmlFor='probe_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Label("Paste up to 6 Probe sequences", htmlFor='probe_input',
                style={"font-size": "15px", "marginLeft": "22px"}),
    html.Div([                                                          # Probe Fasta input box 
     dcc.Textarea(id='probe_input', rows=10, cols=90, value='', required=True,
                placeholder=f'>probe name\nACGTACGT...\n>probe2 name\nACGTACGT...\n',)],
                style={"marginBottom": "20px"}),

    html.Div([                                                          # Restriction enzyme input title
     html.Label("3) Select one or more Restriction enzymes ", htmlFor='rest_enz_input',
                style={"font-size": "20px"}),
     html.Label("(Check out a list of restriction enzymes ", htmlFor='rest_enz_input',
                style={"font-size": "15px"}),
     html.A('Here', href='https://rebase.neb.com/rebase/link_bionetc', target='_blank',
            style={"font-size": "15px", "color": "#007bff", "textDecoration": "none"}),
    html.Label(")", style={"font-size": "15px"})           
                ]),
     html.Div(
     html.Label("Multiple restriction digests are performed sequentially, following the order in which the restriction enzymes are added.", htmlFor='rest_enz_instructions1',
                style={"font-size": "15px", "marginLeft": "22px"})),
    html.Div(id='rest_enz_input',children=[                             # Restriction enzyme input box
         dcc.Input( id= {'type': 'enzyme', 'index': 0}, type='text', value='',
                   placeholder='Enter one restriction enzyme name',
                    style={'width': '400px',"marginBottom": "10px","marginTop": "5px", "display": "block"}),],
                    style={"display": "flex", "flexDirection": "column"}),
    html.Div([                                                          # Add enzyme button  
     html.Button('Add Restriction Enzymes', id='add_enzyme_button', n_clicks=0,
                 style={"font-size": "12px", "marginTop": "10px", 'marginLeft': '250px'})]),
                 
    html.Div(id="enzyme_suggestion_output", style={"color": "orange", "cursor": "pointer"}),      #Output containing enzyme suggestions
    html.Div(id='error_output', style={'color': 'red', 'fontWeight': 'bold','marginTop': '30px'}),      #Output containing the error messages
    
    html.Div([                                                          # Submission button  
     html.Button('Submit', id='submit_button', n_clicks=0,
                 style={"font-size": "20px", "marginTop": "50px", 'marginLeft': '560px'})]),
    
     html.Div([                                         #Gel image
         html.Img(id = 'gel_image', src = '')],
         id='plot_div'),
     html.Div([                                         #HTML
        html.Iframe(id="html-frame", srcDoc= '', style={"width": "100%", "height": "600px"})
     ]),

     html.Div([                                         #Data Storage
         dcc.Store(id = 'store_digested_dict', data = {}),
         dcc.Store(id = 'store_probes_dict', data = {})])
])

########################################################################################################################
########################################################################################################################

# Callback to add more enzyme input boxes for every click on the Add enzyme button

@app.callback(
    Output('rest_enz_input', 'children'),
    Input('add_enzyme_button', 'n_clicks'),
    State('rest_enz_input', 'children')
)

def display_additional_enzyme_inputs(n_clicks, children):
#Each time the Add Enzyme button is clicked, add a new input box.
    if n_clicks == 0:
        return children     # do nothing if not triggered by the button
    if n_clicks < 2:        # limit to max 2 enzyme boxes for now
        new_id = {'type': 'enzyme', 'index': counter['n']}
        new_box = dcc.Input(id= new_id, type='text', value='',
                   style={'width': '400px',"marginBottom": "10px", "display": "block"})
        counter['n'] += 1
        children.append(new_box)
    return children

# CALLBACK: suggests enzyme names as user types AND get the list of the enzyme patterns
@app.callback(
    Output("enzyme_suggestion_output", "children"),
    Input("submit_button", "n_clicks"),
    State({'type': 'enzyme', 'index': ALL}, 'value'),
)
def check_enzyme_suggestions(n_clicks, enzyme_values):
    message = []
    count = 0
    if n_clicks == 0:
        return no_update
    else:
        for msg in re_match(enzyme_values)[1]:
            if '^' in msg:
                message.append(html.Div(
                f"{enzyme_values[count]}: {msg}",
                style={"marginBottom": "8px", 'color': 'green'}
                ))
            else:
                message.append(html.Div(
                f"{enzyme_values[count]}: {msg}",
                style={"marginBottom": "8px"}
                ))
            count += 1
        return message

########################################################################################################################
########################################################################################################################


# CALLBACK: Capture user inputs and store in variables when submit button is clicked

@app.callback(
    Output('store_digested_dict', 'data'), #this will be enzyme_values
    Output('store_probes_dict', 'data'),#this will be probe_values
    Output('error_output', 'children'),
    Input('submit_button', 'n_clicks'), # triggers when Submit button is clicked
    State('Upload button', 'contents'),
    State('fasta_input', 'value'),  # reads FASTA input
    State({'type': 'enzyme', 'index': ALL}, 'value'), # read values of all enzyme boxes
    State('probe_input', 'value')  # reads Probe FASTA input
)

def capture_inputs(n_clicks, fasta_values_upload, fasta_values_text, enzyme_values, probe_values):
    if n_clicks == 0:
        return no_update, no_update, no_update # Do nothing if not triggered by the button
    fasta_values = ''
    error_msg = []
###########################################

    if not fasta_values_upload and not fasta_values_text:
        error_msg.append("Error:Invalid FASTA format for DNA sequence. Please ensure your input starts with '>'.")
    if fasta_values_upload and fasta_values_text:
        error_msg.append("Error: Please EITHER upload OR paste your FASTA file.")
    if fasta_values_upload:
        content_type, content_string = fasta_values_upload.split(',')
        decoded_bytes = base64.b64decode(content_string)
        text = decoded_bytes.decode('utf-8')
        text = text.lstrip('\ufeff')
        fasta_values = text
        #print(f"x{fasta_values[:100]}")
    if fasta_values_text:
        fasta_values = fasta_values_text
        #print(f"y{fasta_values[:100]}")
 #########################################   
    if fasta_values == '':
        error_list = []
        for msg in error_msg:
            error_list.append(html.Div(
            f"{msg}",
            style={"marginBottom": "8px", 'color': 'red'}
            ))
        return no_update, no_update, error_list
    else:
        if not fasta_values.startswith('>'): #check the FASTA has header
            error_msg.append("Error: Invalid FASTA format for DNA sequence. Please ensure your input starts with '>'.")
        for line in fasta_values.split('\n'): #Check only valid characters in FASTA sequence
            if '>' not in line:  # Skip header lines
                line = line.upper().strip()
                if not set(line).issubset(['A','C','G','T','N']):
                    error_msg.append("Error: Invalid characters in DNA FASTA sequence. Only A, C, G, T, N are allowed.")
        
        if not probe_values or not probe_values.startswith('>'):  #check at least one probe entered
            error_msg.append("Error: Invalid FASTA format for Probe sequence. Please ensure your input starts with '>'.")
        for line in probe_values.split('\n'): #Check only valid characters in FASTA sequence
            if not line.startswith('>'):  # Skip header lines
                line = line.upper().strip()
                if not set(line).issubset(['A','C','G','T','N','R','Y']):
                    error_msg.append("Error: Invalid characters in Probe FASTA sequence. Only A, C, G, T, N, R, Y are allowed.")
        
        if not enzyme_values or all(e is None or e.strip() == "" for e in enzyme_values):  #check at least one enzyme entered. 
            error_msg.append("Error: Please enter one restriction enzyme.")
        for enz in enzyme_values:
            if ' ' in enz:
                error_msg.append("Error: Enzyme names cannot contain spaces.")
        
        if error_msg:
            error_list = []
            for msg in error_msg:
                error_list.append(html.Div(
                f"{msg}",
                style={"marginBottom": "8px", 'color': 'red'}
                ))
            return no_update, no_update, error_list
        
        if len(re_match(enzyme_values)[0]) == len(enzyme_values):  #check all enzymes valid
            if not error_msg:  # proceed only if no errors so far
                fasta_dict = read_fasta(fasta_values) #fasta dictionary
                #(fasta_dict)
                probes_dict = read_fasta(probe_values)
                motifs = re_match(enzyme_values)[0] #enzyme patterns list
                digested_dict = re_digest(fasta_dict, motifs)
                #print(digested_dict)

                return digested_dict, probes_dict, []

########################################################################################################################
########################################################################################################################

#Generating gel PNG image
@app.callback(
Output("gel_image", 'src'),
Input('store_digested_dict', 'data'),
State('fasta_input', 'value'),
State('Upload button', 'contents'),
State({'type': 'enzyme', 'index': ALL}, 'value'),
State('probe_input', 'value')
)
def update_gel_image(n_clicks, fasta_values_text, fasta_values_upload, enzyme_values, probe_values):
    #Graphical parameters
    ladder = {'25bp':25 , '50bp':50, '100bp': 100, '200bp': 200, '500bp': 500,
        '1000bp': 1000, '2000bp': 2000, '3000bp': 3000, '5000bp' : 5000 , '10000bp':10000 , '20000bp':20000 , '30000bp':30000 , '100000bp':100000}
    probe_colors = ["yellow" , "lime" , "magenta" , "cyan" , "orange", "teal" , "blue", "coral" , "darkgreen" , "tan", "sienna", "plum", "lavendar"]

    if n_clicks == 0:
        return no_update
    
    if fasta_values_upload:
            content_type, content_string = fasta_values_upload.split(',')
            decoded_bytes = base64.b64decode(content_string)
            text = decoded_bytes.decode('utf-8')
            text = text.lstrip('\ufeff')
            fasta_values = text
            #print(f"x{fasta_values[:100]}")
    if fasta_values_text:
        fasta_values = fasta_values_text

    # Process inputs as before
    fasta_dict = read_fasta(fasta_values)
    probes_dict = read_fasta(probe_values)
    motifs = re_match(enzyme_values)[0]
    digested_dict = re_digest(fasta_dict, motifs)

    #Create URI image
    gel_figure = gel_image(ladder, digested_dict, probes_dict, probe_colors)
    return gel_figure

########################################################################################################################
########################################################################################################################

#Generating HTML output html-frame
@app.callback(
Output("html-frame", 'srcDoc'),
Input('store_digested_dict', 'data'),
State('fasta_input', 'value'),
State('Upload button', 'contents'),
State({'type': 'enzyme', 'index': ALL}, 'value'),
State('probe_input', 'value')
)
def update_html(n_clicks, fasta_values_text, fasta_values_upload, enzyme_values, probe_values):

    probe_colors = ["yellow" , "lime" , "magenta" , "cyan" , "orange", "teal" , "blue", "coral" , "darkgreen" , "tan", "sienna", "plum", "lavendar"]

    if n_clicks == 0:
        return no_update
    
    if fasta_values_upload:
            content_type, content_string = fasta_values_upload.split(',')
            decoded_bytes = base64.b64decode(content_string)
            text = decoded_bytes.decode('utf-8')
            text = text.lstrip('\ufeff')
            fasta_values = text
            #print(f"x{fasta_values[:100]}")
    if fasta_values_text:
        fasta_values = fasta_values_text

    # Process inputs as before
    fasta_dict = read_fasta(fasta_values)
    probes_dict = read_fasta(probe_values)
    motifs = re_match(enzyme_values)[0]
    digested_dict = re_digest(fasta_dict, motifs)

    #html_display = generate_html(digested_dict, probes_dict, probe_colors)
    html_content = generate_html(digested_dict, probes_dict, probe_colors)
    #html_raw = base64.b64decode(html_content.split(",")[1]).decode('utf-8')

    i = '\n'.join(html_content)
    return i
    
   

########################################################################################################################
########################################################################################################################

# #Generating txt file
# save_fragment_lengths("TEXT_output.txt", digested_dict)





if __name__ == "__main__":
    app.run(debug=True)