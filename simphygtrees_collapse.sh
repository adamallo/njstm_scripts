#!/bin/bash
declare -a positionals
while [[ $# > 0 ]]
do
	key="$1"

	case $key in
    		-r|--remove)
    		remove=1
    		shift # past argument
   		;;
    		*)
		positionals+=($1)
		shift
    		;;
	esac
done

if [[ ${#positionals[@]} -ne 1 ]] #Only one positional
then
	usage=1
else
	dir=${positionals[0]}
	if [[ ! -d $dir ]]
	then
		usage=1
	else
		usage=0
	fi
fi
if [[ $usage -eq 1 ]]
then
	echo "Usage $0 SimPhy_outputfolder [-r|--remove]"
	echo "The remove option removes the original input files"
	exit
fi

for i in $dir/*
do 
	if [[ -d $i ]]
	then 
		cat $i/g_trees* >> $i/original_g_trees.tree
		if [[ $remove -eq 1 ]]
		then
			rm -f $i/g_trees*
		fi
	fi
done
