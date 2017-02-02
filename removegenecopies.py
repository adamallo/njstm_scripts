#!/opt/local/bin/python

### Imports ###
import dendropy
from dendropy import TreeList,Tree
import sys
import argparse
from os import walk
import glob
#from scipy.stats import truncnorm
import numpy
from numpy.random import normal
from joblib import Parallel, delayed

### My functions ###
def remove_taxa_prov(r,nodelist,prov_removal):
	discarded=list()
	rs=r.random_sample(len(nodelist))
	for i in xrange(len(nodelist)):
		if rs[i]<= prov_removal:
			 discarded.append(nodelist[i].taxon) 
	return discarded

def remove_taxa_tagprobs(r,nodelist,prob_dict):
	discarded=list()
	rs=r.random_sample(len(nodelist))
	for i in xrange(len(nodelist)):
		assert prob_dict[nodelist[i].taxon.label], "Pobability not found!\n"
		if rs[i]<= prob_dict[nodelist[i].taxon.label]:
			discarded.append(nodelist[i].taxon)
	return discarded

def truncated_normal(r,n,mean,sd,min=None,max=None):
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
    			rnum=r.normal(mean,sd)
    			if (min is None) or (min <= rnum):
      				accept+=1
    			if (max is None) or (max >= rnum):
      				accept+=1
    			if accept == 2:
      				rnumbers.append(rnum)
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
parser.add_argument("-ncores",type=int,help="Number of cores")
args = parser.parse_args()

###Random number machinery initialization
if args.s:
	seed=args.s
else:
	seed=numpy.random.randint()

numpy.random.seed(seed)
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

if len(folders) == 0:
	raise NameError("0 Detected folders")

###For each replicate
def main (folder=None,seed=None):
	print("Folder %s, seed %s") % (folder,seed)
	r=numpy.random.RandomState(seed)
	gene_trees=TreeList()
	taxa = dendropy.TaxonNamespace()
	treefiles=glob.glob(args.sd+"/"+folder+"/g_trees*.trees")
	tree_yielder=Tree.yield_from_files(files=treefiles,schema="newick",rooting="default-rooted",preserve_underscores=True,taxon_namespace=taxa)
	#Modify gene trees
	#I have to modify here the trees
	if args.mk=="random":
		for gtree in tree_yielder:
			onodes=gtree.leaf_nodes()
			nodes=remove_taxa_prov(r,onodes,args.pr)
			if len(nodes) < len(onodes)-3: #Tree with missing leaves
				gtree.prune_taxa(nodes,update_bipartitions=False, suppress_unifurcations=True)
				gene_trees.append(gtree)
			else:	#The whole tree is missing (the tree would have 3 leaves or less, which is not an unrooted tree)
				continue
	elif args.mk=="byindividual":
		tagProbs=None
		for gtree in tree_yielder:
                        onodes=gtree.leaf_nodes()
			if not tagProbs:
				tagProbs={}
				probs=truncated_normal(r,n=len(onodes),mean=args.pr,sd=args.ist,min=args.itmin,max=args.itmax) #one prob for each leaf
				for leafi in xrange(len(onodes)):
					tagProbs[onodes[leafi].taxon.label]=probs[leafi]#assigment to leaf labels in the dictionary
                        nodes=remove_taxa_tagprobs(r,onodes,tagProbs)
			if len(nodes) < len(onodes)-3: #Tree with missing leaves
                                gtree.prune_taxa(nodes,update_bipartitions=False, suppress_unifurcations=True)
                                gene_trees.append(gtree)
                        else:   #The whole tree is missing (the tree would have 3 leaves or less, which is not an unrooted tree)
                                continue
	else:
		print("Yet unsupported option")
	#Write gene trees
	gene_trees.write(path=args.sd+"/"+folder+"/"+args.o,schema="newick")

seeds= numpy.random.randint(numpy.iinfo(numpy.int32).max,size=len(folders))

if args.ncores > 1:
	Parallel(n_jobs=args.ncores) (delayed(main)(folder=folders[i],seed=seeds[i]) for i in xrange(len(folders)))
else:
	[main(folder=folders[i],seed=seeds[i]) for i in xrange(len(folders))]

print("Done!")
