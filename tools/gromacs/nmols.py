#!/usr/bin/env python3

'#extracts the number of molecules succesfully added by gmx insert-molecules.'

import re


inFile = open('verbose.txt', 'r')
outFile = open("addedmols.txt", "w")
lines = inFile.read()
result = re.search('Added(.*)molecules', lines)

inFile.close()
outFile.write(str(result.group(1)))
outFile.close()
