#!/bin/bash

#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

declare -a positionals
id="*"
while [[ $# > 0 ]]
do
	key="$1"

	case $key in
    		-r|--remove)
    		remove=1
    		shift # past argument
   		;;
		-i|--id)
		id=$2
		shift 2
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
	echo "Usage $0 SimPhy_outputfolder [-r|--remove] [-i|--id id]"
	echo "The remove option removes the original input files"
	exit
fi

if [ "$SGE_TASK_ID" != "" ] && [ "$SGE_TASK_ID" != "undefined" ]
then
        id=$(printf "%06d" $SGE_TASK_ID)
fi

cd $dir

for i in $id
do 
	if [[ -d $i ]]
	then 
		echo "working on $i"
		cat $i/g_trees* >> $i/original_g_trees.tree
		if [[ $remove -eq 1 ]]
		then
			rm -f $i/g_trees*
		fi
	fi
done
