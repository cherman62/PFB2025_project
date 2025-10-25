from dash import Dash, html, dcc, Input, Output, State, ALL, no_update
import pandas as pd
from read_fasta import read_fasta


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
    html.Div([                                                          # Fasta input title  
     html.Label("Paste one or more FASTA sequences", htmlFor='fasta_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Div([                                                          # Fasta input box 
     dcc.Textarea(id='fasta_input', rows=20, cols=90, required=True,
                placeholder=f'>sequence name\nACGTACGT...\n>sequence2 name\nACGTACGT...\n',)],
                style={"marginBottom": "20px"}),
    html.Div([                                                          # Restriction enzyme input title
     html.Label("Restriction enzyme", htmlFor='rest_enz_input',
                style={"font-size": "20px", "marginRight": "10px"})]),
    html.Div(id='rest_enz_input',children=[                             # Restriction enzyme input box
         dcc.Input( id= {'type': 'enzyme', 'index': 0}, type='text', required= True,
                    style={'width': '250px',"marginBottom": "10px", "display": "block"}),],
                    style={"display": "flex", "flexDirection": "column"}),
    html.Div([                                                          # Add enzyme button  
     html.Button('Add Restriction Enzymes', id='add_enzyme_button', n_clicks=0,
                 style={"font-size": "10px", "marginTop": "10px", 'marginLeft': '200px'})]),
    html.Div([                                                          # Add enzyme button  
     html.Button('Submit', id='submit_button', n_clicks=0,
                 style={"font-size": "10px", "marginTop": "100px", 'marginLeft': '600px'})]),
    html.Div(id="seq_output"),
    html.Div(id="enz_output")
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
    
    new_id = {'type': 'enzyme', 'index': counter['n']}
    new_box = dcc.Input(id= new_id, type='text',
                   style={'width': '250px',"marginBottom": "10px", "display": "block"})
    counter['n'] += 1
    children.append(new_box)
    return children


# CALLBACK: Capture user inputs and store in variables

@app.callback(
    Output('seq_output', 'children'), #since this is the first output it will be fasta_input
    Output('enz_output', 'children'), #this will be enzyme_values
    Input('submit_button', 'n_clicks'), # triggers when Submit button is clicked
    State('fasta_input', 'value'),  # reads FASTA input
    State({'type': 'enzyme', 'index': ALL}, 'value') # read values of all enzyme boxes
)

def capture_inputs(n_clicks, fasta_input, enzyme_values):
    if n_clicks == 0:
        return no_update, no_update  # Do nothing if not triggered by the button
    if '>' not in fasta_input: #check the FASTA has header
        return "Error: Invalid FASTA format. Please ensure your input starts with '>'.", no_update
    for line in fasta_input.split('\n'): #Check only valid characters in FASTA sequence
        print(line)
        if '>' not in line:  # Skip header lines
            line = line.upper().strip()
            if not set(line).issubset(['A','C','G','T','N']):
                return "Error: Invalid characters in FASTA sequence. Only A, C, G, T, N are allowed.", no_update
    if enzyme_values == [None]:  #check at least one enzyme entered
        return "Error: Please enter one restriction enzyme.", no_update
    for enz in enzyme_values:
        if ' ' in enz:
            return "Error: Enzyme names cannot contain spaces.", no_update
    else:
        read_fasta(fasta_input)

if __name__ == "__main__":
    app.run(debug=True)