#!/opt/local/bin/python

### Imports ###
import dendropy
from dendropy import TreeList,Tree
import sys
import argparse
from os import walk
import glob
import random
#from scipy.stats import truncnorm
from numpy.random import normal
from joblib import Parallel, delayed

### My functions ###
def remove_taxa_prov(nodelist,prov_removal):
	discarded=list()
	for node in nodelist:
		if random.random()<= prov_removal:
			 discarded.append(node.taxon) 
	return discarded

def remove_taxa_tagprobs(nodelist,prob_dict):
	discarded=list()
	for node in nodelist:
		assert prob_dict[node.taxon.label], "Pobability not found!\n"
		if random.random()<= prob_dict[node.taxon.label]:
			discarded.append(node.taxon)
	return discarded

def truncated_normal(n,mean,sd,min=None,max=None):
#	if (min is not None) and (max is not None):
#		return [truncnorm((min-mean)/sd, (max-mean)/sd) for i in xrange(n)]
#	else:
		maxit=1000000
		rnumbers=list()
  		i=0
		it=0
  		while (i<n):
			assert it<=maxit, "Maximum iteration sampling the truncated normal reached\n"
    			accept=0
    			r=normal(mean,sd)
    			if (min is None) or (min <= r):
      				accept+=1
    			if (max is None) or (max >= r):
      				accept+=1
    			if accept == 2:
      				rnumbers.append(r)
      				i+=1
			it+=1
		return rnumbers


### Main ###

### Argparse
parser = argparse.ArgumentParser(description="Generates missing data removing gene copies from gene trees simulated using SimPhy",prog="removegenecopies.py")
parser.add_argument("sd",type=str,help="Main simphy output directory",metavar="inputdir")
parser.add_argument("pr",type=float,help="Missing data probability",metavar="probability")
parser.add_argument("o",type=str,help="Output file name",metavar="outname")
parser.add_argument("-mk",type=str,default="random",help="Missing data generation scheme. Implemented schemes: random and byindividual",choices=["random","byindividual","bygene"],metavar="data_scheme")
parser.add_argument("-ist",type=float,help="Standard deviation for the truncated normal distribution with mean p to sample by-individual missing probabilities, only used if -mk byindividual",metavar="standardeviation",default=1)
parser.add_argument("-itmin",type=float,help="Minimum value to truncate the normal distribution to sample by-individual missing probabilities, only used if -mk byindividual",metavar="truncmin",default=None)
parser.add_argument("-itmax",type=float,help="Maximum value to truncate the normal distribution to sample by-individual missing probabilities,, only used if -mk byindividual",metavar="truncmax",default=None)
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
	elif args.mk=="byindividual":
		tagProbs=None
		for gtree in tree_yielder:
                        onodes=gtree.leaf_nodes()
			if not tagProbs:
				tagProbs={}
				probs=truncated_normal(n=len(onodes),mean=args.pr,sd=args.ist,min=args.itmin,max=args.itmax) #one prob for each leaf
				for leafi in xrange(len(onodes)):
					tagProbs[onodes[leafi].taxon.label]=probs[leafi]#assigment to leaf labels in the dictionary
                        nodes=remove_taxa_tagprobs(onodes,tagProbs)
			if len(nodes) < len(onodes)-1: #Tree with missing leaves
                                gtree.prune_taxa(nodes,update_bipartitions=False, suppress_unifurcations=True)
                                gene_trees.append(gtree)
                        else:   #The whole tree is missing
                                continue
	else:
		print("Yet unsupported option")
	#Write gene trees
	gene_trees.write(path=args.sd+"/"+folder+"/"+args.o,schema="newick")
print("Done!")
