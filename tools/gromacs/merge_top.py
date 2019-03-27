import sys
import re

def combine_tops(top_text, itp_texts):
    """
    Search through parent topology top_text and replace #include lines with the relevant child topologies from the dictionary itp_texts
    """
    for itp in itp_texts:
        spl = re.split('#include ".*{}"\n'.format(itp), top_text) # split on include string
        top_text = itp_texts[itp].join(spl)
    return top_text

top = sys.argv[1] # parent topology file
itps_file = sys.argv[2] # file with list of child topologies (.itp files)

with open(itps_file) as f:
    itps = f.read().split()

with open(top, 'r') as f:
    top_text = f.read()

itp_texts = {} # create dictionary of child topologies
for itp in itps:
    with open(itp, 'r') as f:
        itp_texts[itp] = f.read()

for itp in itp_texts:
    # child tops may also refer to each other; we need to check this
    itp_texts[itp] = combine_tops(itp_texts[itp], itp_texts)

with open('top_output.top', 'w') as f:
    f.write(combine_tops(top_text, itp_texts)) # now combine all children into the parent