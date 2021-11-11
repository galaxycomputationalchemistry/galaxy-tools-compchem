#!/usr/bin/env python3

'#this script extracts the number of molecules succesfully added by gmx insert-molecules, and saves it as a string in a text file. This value can then be used to update topology information later on before running a simulation.'

import re


inFile = open('verbose.txt', 'r')
outFile = open("addedmols.txt", "w")
lines = inFile.read()
result = re.search('Added(.*)molecules', lines)

inFile.close()
outFile.write(str(result.group(1)))
outFile.close()
