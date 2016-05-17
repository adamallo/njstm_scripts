#!/opt/local/bin/python

### Imports ###
import dendropy
from dendropy import TreeList,Tree
import sys
import argparse
from os import walk
import glob
import random

### My functions ###
def remove_taxa_prov(nodelist,prov_removal):
	discarded=list()
	for node in nodelist:
		if random.random()<= prov_removal:
			 discarded.append(node.taxon) 
	return discarded
### Main ###

### Argparse
parser = argparse.ArgumentParser(description="Generates missing data removing gene copies from gene trees simulated using SimPhy",prog="removegenecopies.py")
parser.add_argument("sd",type=str,help="Main simphy output directory",metavar="inputdir")
parser.add_argument("pr",type=float,help="Missing data probability",metavar="probability")
parser.add_argument("o",type=str,help="Output file name",metavar="outname")
parser.add_argument("-mk",type=str,default="random",help="Missing data generation scheme",choices=["random","2sp","bygene"],metavar="data_scheme")
parser.add_argument("-a",type=float,help="Alpha parameter for the by-gene gamma distribution, only used if -mk bygene",metavar="alpha",default=1)
parser.add_argument("-s",type=int,help="Random number generator seed",metavar="seed")
args = parser.parse_args()

###Random number machinery initialization
if args.s:
	seed=args.s
else:
	seed=random.randint(0,sys.maxint)

random.seed(seed)
print("Seed: %d" % seed)

###Some variables to recycle
nodes=list()
onodes=list()
id=0
###Get replicate folders
folders=list()
for (dirpath, dirnames, filenames) in walk(args.sd):
        folders.extend(dirnames)
        break

###For each replicate
for folder in folders:
	#Get gene trees	
#	for treefile in glob.glob(args.sd+"/"+folder+"/g_trees*.trees"):
#		gene_trees.append(Tree.get(path=treefile,schema="newick",rooting="default-rooted",preserve_underscores=True))
	gene_trees=TreeList()
	taxa = dendropy.TaxonNamespace()
	treefiles=glob.glob(args.sd+"/"+folder+"/g_trees*.trees")
	tree_yielder=Tree.yield_from_files(files=treefiles,schema="newick",rooting="default-rooted",preserve_underscores=True,taxon_namespace=taxa)
	#Modify gene trees
	#I have to modify here the trees
	if args.mk=="random":
		for gtree in tree_yielder:
			onodes=gtree.leaf_nodes()
			nodes=remove_taxa_prov(onodes,args.pr)
			if len(nodes) < len(onodes)-1: #Tree with missing leaves
				gtree.prune_taxa(nodes,update_bipartitions=False, suppress_unifurcations=True)
				gene_trees.append(gtree)
			else:	#The whole tree is missing
				continue
	else:
		print("Yet unsupported option")
	#Write gene trees
	gene_trees.write(path=args.sd+"/"+folder+"/"+args.o,schema="newick")
print("Done!")
