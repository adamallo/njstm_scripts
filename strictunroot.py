#!/opt/local/bin/python

### Imports ###
import dendropy
from dendropy import TreeList,Tree
import sys
import argparse
from os import walk
import glob


### Main ###

### Argparse
parser = argparse.ArgumentParser(description="Reads a newick trees and reroots it with a basal trifurcation",prog="strictunroot.py")
parser.add_argument("-i",required=True,type=str,help="Input newick tree name")
parser.add_argument("-o",required=True,type=str,help="Output file name")
args = parser.parse_args()

###Main
tree_yielder=Tree.yield_from_files(files=[args.i],schema="newick",rooting="default-rooted",preserve_underscores=True)
for tree in tree_yielder:
    tree.collapse_basal_bifurcation()
    tree.write(path=args.o,schema="newick",unquoted_underscores=True,suppress_rooting=True)
print("Done!")
