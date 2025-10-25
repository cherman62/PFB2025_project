#!/usr/bin/env python

from dash import Dash, html, dcc, Input, Output, State, ALL, no_update
import re
from read_fasta import read_fasta
from read_re import re_match


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
     html.Div([                                                          # Header
    html.Header("Some kind of description here.................", 
                 style={"font-size": "20px"})],
     style={"marginBottom": "30px"}),
    html.Div([                                                          # Fasta input title  
     html.Label("1) Paste one or more DNA sequences in FASTA format", htmlFor='fasta_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Div([                                                          # Fasta input box 
     dcc.Textarea(id='fasta_input', rows=20, cols=90, required=True,
                placeholder=f'>sequence name\nACGTACGT...\n>sequence2 name\nACGTACGT...\n',)],
                style={"marginBottom": "20px"}),
                html.Div([                                              # Probe Fasta input title  
     html.Label("2) Paste one or more Probe sequences in FASTA format", htmlFor='probe_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Div([                                                          # Probe Fasta input box 
     dcc.Textarea(id='probe_input', rows=10, cols=90, required=True,
                placeholder=f'>probe name\nACGTACGT...\n>probe2 name\nACGTACGT...\n',)],
                style={"marginBottom": "20px"}),
    html.Div(                                                          # Restriction enzyme input title
     html.Label("3) Select one or more Restriction enzymes", htmlFor='rest_enz_input',
                style={"font-size": "20px", "marginRight": "10px"})),
    html.Div(id='rest_enz_input',children=[                             # Restriction enzyme input box
         dcc.Input( id= {'type': 'enzyme', 'index': 0}, type='text', required= True,
                    style={'width': '400px',"marginBottom": "10px", "display": "block"}),],
                    style={"display": "flex", "flexDirection": "column"}),
    html.Div([                                                          # Add enzyme button  
     html.Button('Add Restriction Enzymes', id='add_enzyme_button', n_clicks=0,
                 style={"font-size": "12px", "marginTop": "10px", 'marginLeft': '250px'})]),
    html.Div(id="enzyme_suggestion_output", style={"color": "orange", "cursor": "pointer"}),      #Output containing enzyme suggestions
    html.Div(id='error_output', style={'color': 'red', 'fontWeight': 'bold','marginTop': '30px'}),      #Output containing the error messages
    html.Div([                                                          # Submission button  
     html.Button('Submit', id='submit_button', n_clicks=0,
                 style={"font-size": "20px", "marginTop": "50px", 'marginLeft': '560px'})]),
    
    html.Div(id="seq_output"),
    html.Div(id="enz_output"),
    html.Div(id="probe_output")
])


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
        new_box = dcc.Input(id= new_id, type='text',
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
    motif_list = []
    message = []
    if n_clicks == 0:
        return no_update
    else:
        for enz in enzyme_values:
            if enz:
                suggestion = re_match(enz)
                motif_list.append(re_match(enz))
                if '^' in suggestion:
                    message.append(html.Div(
                    f"{enz}: {suggestion}",
                    style={"marginBottom": "8px", 'color': 'green'}
                ))
                else:
                    message.append(html.Div(
                    f"{enz}: {suggestion}",
                    style={"marginBottom": "8px"}
                ))
        return message

# CALLBACK: Capture user inputs and store in variables

@app.callback(
    Output('error_output', 'children'),  # Output for error messages
    Output('seq_output', 'children'), #since this is the first output it will be fasta_values
    Output('enz_output', 'children'), #this will be enzyme_values
    Output('probe_output', 'children'),#this will be probe_values
    Input('submit_button', 'n_clicks'), # triggers when Submit button is clicked
    State('fasta_input', 'value'),  # reads FASTA input
    State({'type': 'enzyme', 'index': ALL}, 'value'), # read values of all enzyme boxes
    State('probe_input', 'value')  # reads Probe FASTA input
)

def capture_inputs(n_clicks, fasta_values, enzyme_values, probe_values):
    if n_clicks == 0:
        return no_update, no_update, no_update, no_update  # Do nothing if not triggered by the button
    if not fasta_values.startswith('>'): #check the FASTA has header
        return ("Error: Invalid FASTA format for DNA sequence. Please ensure your input starts with '>'.", no_update, no_update, no_update)
    for line in fasta_values.split('\n'): #Check only valid characters in FASTA sequence
        if '>' not in line:  # Skip header lines
            line = line.upper().strip()
            if not set(line).issubset(['A','C','G','T','N']):
                return ("Error: Invalid characters in DNA FASTA sequence. Only A, C, G, T, N are allowed.", no_update, no_update, no_update)
    if enzyme_values == [None]:  #check at least one enzyme entered. !!!!!NOT WORKING YET!!!!!
        return ("Error: Please enter one restriction enzyme.", no_update, no_update, no_update)
    for enz in enzyme_values:
        if ' ' in enz:
            return ("Error: Enzyme names cannot contain spaces.", no_update, no_update, no_update)
    if not probe_values.startswith('>'):  #check at least one probe entered
        return ("Error: Invalid FASTA format for Probe sequence. Please ensure your input starts with '>'.", no_update, no_update, no_update)
    for line in probe_values.split('\n'): #Check only valid characters in FASTA sequence
        if '>' not in line:  # Skip header lines
            line = line.upper().strip()
            if not set(line).issubset(['A','C','G','T','N','R','Y']):
                return ("Error: Invalid characters in Probe FASTA sequence. Only A, C, G, T, N, R, Y are allowed.", no_update, no_update, no_update)
    else:
        read_fasta(fasta_values)
        read_fasta(probe_values)
        print(type(enzyme_values))
        return (None, fasta_values, enzyme_values, probe_values) #Temporaryly return the inputs for testing
        

if __name__ == "__main__":
    app.run(debug=True)