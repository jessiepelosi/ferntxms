'''
Purpose: Summarize results from MAPS simulations for runFisher_null.pl and runFisher_positive.pl
Usage: python summarize_MAPS.py [number of nodes] [Number of trees analyzed in each MAPS analysis] [list name]

This should be run in the same directory as the CSV outputs from MAPS for the null and positive simulations. Only run
this on the null or positive simulations independently (one directory for null and one directory for positive). The output
from this script can be used as the input for the runFisher_null.pl and runFisher_positive.pl scripts (or modifications
of those).

Jessie Pelosi, last modified Nov 11, 2021
'''

import pandas as pd
import glob
import sys
from pprint import pprint
args = sys.argv
NoNode = sys.argv[1]
NoTrees = sys.argv[2]
ListName = sys.argv[3]

sim_files = glob.glob("./*subtree.csv")

# Generate an empty list of dictionaries

dictionaries = {}

print("Generating Dictionaries.")

# Create nested dictionaries with four metrics for each node

for Node in range(int(NoNode)):
    NodeVal = 'N%i' % Node
    dictionaries[NodeVal]={"NonDupPer": [], "DupPer": [], "NonDupNo": [], "DupNo": []}

#pprint(dictionaries) # Test

print("Assembling Data.")

for file in sim_files:
        df = pd.read_csv(file, index_col=None, header=0)
        for i in range(int(NoNode)):
                #print(i)
                #print(type(i))
                dictionaries[("N%i" % i)]['NonDupPer'].append(float(df.iat[i,1].rstrip("%"))) #Remove % sign here for downstream stuff
                dictionaries[("N%i" % i)]['DupPer'].append(float(df.iat[i,2].rstrip("%"))) #Remove % sign here for downstream stuff
                dictionaries[("N%i" % i)]['NonDupNo'].append(float(df.iat[i,3]))
                dictionaries[("N%i" % i)]['DupNo'].append(float(df.iat[i,4]))

# Make empty lists for each metric of interest

MRCA = [] # Node
NonDupPer = [] # Percent of non-duplicated subtrees
DupPer = [] # Percent of duplicated subtrees
NonDupNo = [] # Number of non-duplicated subtrees
DupNo = [] # Number of duplicated subtrees

# Get averages for each metric and put into the empty lists for each node

for i in range(int(NoNode)):
        MRCA.append('N%i' % (i + 1))
        NonDupPer.append(sum(dictionaries[("N%i" % i)]['NonDupPer'])/len(dictionaries[("N%i" % i)]['NonDupPer']))
        DupPer.append(sum(dictionaries[("N%i" % i)]['DupPer'])/len(dictionaries[("N%i" % i)]['DupPer']))
        NonDupNo.append(sum(dictionaries[("N%i" % i)]['NonDupNo'])/len(dictionaries[("N%i" % i)]['NonDupNo']))
        DupNo.append(sum(dictionaries[("N%i" % i)]['DupNo'])/len(dictionaries[("N%i" % i)]['DupNo']))

#print(NonDupPer)
#print(DupPer)
#print(NonDupNo)
#print(DupNo)

# Create dataframe from lists, edit, and write out to csv to use in runFisher*.pl scripts

temp_df = pd.DataFrame((list(zip(MRCA,NonDupPer, DupPer, NonDupNo, DupNo))), columns = ['MRCA','NonDupPer', 'DupPer', 'NonDupNo', 'DupNo'])
temp_df['NonDupPer'] = temp_df['NonDupPer'].astype(str) +'%' # Add % back, needs to be here for runFisher*.pl scripts
temp_df['DupPer'] = temp_df['DupPer'].astype(str) +"%" # Add % back, needs to be here for runFisher*.pl scripts
temp_df['TotalNoSubtrees'] = temp_df['DupNo']+temp_df['NonDupNo'] # Make new column for total number of subtrees analyzed
temp_df['NoTrees'] = NoTrees # Make new column for total number of trees analyze in each maps analysis
print("Writing CSV.")
temp_df.to_csv(ListName + ".meanMapsOut.csv",index=False) # Write to CSV!

print("Done! Program finished successfully.")
