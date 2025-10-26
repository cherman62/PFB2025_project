#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math

#to visualize fragmented DNA in a gel, at first I don't need to worry about sequences, instead I just want to see where the fragments sizes are on the gel/membrane#
#as an example, I want to define the fragment sizes I was working with in my test file for example: 
fragment_sizes = {"Fragment_1_250bp" : 250 , 
                  "Fragment_2_500bp" : 500 , 
                  "Fragment_3_800bp" : 800 , 
                  "Fragment_4_1200bp" : 1200 , 
                  "Fragment_5_2000bp" : 2000 , 
                  "Fragment_6_3000bp" : 3000}

#DNA ladder(standard sizes in bp)
ladder = {'100bp' : 100,
          '200bp' : 200,
          '500bp' : 500, 
          '1000bp': 1000, 
          '2000bp' : 2000, 
          '3000bp' : 3000}


#now I want to simulate migration distance (i.e., smaller DNA fragments move further)
#here, it is suggested to use an inverse log relationship to approximate real gels 
#why log? in real agarose gels, migration distance is roughly inverstly proportional to the log of DNA size: smaller gragments move farther than the larger fragments. Using a log give a realistic-looking gel spacing. 

def migration_distance(bp, gel_height=10, max_bp=5000, min_bp=100):#here the max_distnace=10 is an optional parameter representing the maximum vertical distance in the plot (i.e., the gel height)
    if bp <=0:
        bp = 1 #logarithmic scaling with a min and max reference 
    return gel_height * (math.log(max_bp) - math.log(bp)) / (math.log(max_bp) - math.log(min_bp))

                #math.log(bp) computes the natural logarithm of the fragment size 
                #dividing by math.log(3000) normalizes the values so that the largest fragment(3000) becomes 1
    #the function will return a y-coordinate for where to place the fragment of the gel. 

#create plot
fig, ax = plt.subplots(figsize=(5, 8))
ax.set_xlim(0, 3)
ax.set_ylim(0,10)
ax.set_yticks([])
ax.set_title('In Silico DNA Gel', color = 'white', fontsize=12, pad=20)
ax.set_facecolor ('black') # set gel background color to black 
fig.patch.set_facecolor('black') #set figure background to black 

#reverse y-axis so largest fragments are at the top
ax.invert_yaxis()

#draw ladder bands in red 
for i, (name, bp) in enumerate(ladder.items()):
    y=migration_distance(bp)
    ax.hlines(y,0.8, 1.2, color = 'red', linewidth =4)
    ax.text(0.1, y, f'{bp}bp', verticalalignment = 'center', color = 'red')

#draw sample bands in cyan 
    for i, (name, bp) in enumerate(fragment_sizes.items()):
        y = migration_distance(bp)
        ax.hlines(y, 1.8, 2.2, color = 'cyan' , linewidth=4)
        ax.text(2.4, y, name, verticalalignment = 'center', color = 'cyan')

#add lane labels at the top of the gel 
ax.text(1.0,-0.03,"DNA Ladder", color ='red', horizontalalignment='center' , fontsize=12, fontweight = 'bold')
ax.text(2.0,-0.03,"Sample", color ='cyan', horizontalalignment='center' , fontsize=12, fontweight = 'bold')


#display gel 
plt.savefig("in_silico_test_gel.png" , dpi = 300, bbox_inches='tight', facecolor = fig.get_facecolor())
plt.show()

