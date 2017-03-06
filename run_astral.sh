#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH

H=/home/dmallo/njstM/data/broadsims/
nDigits=5

if [ "$SGE_TASK_ID" == "" ] || [ "$SGE_TASK_ID" == "undefined" ]
then
	id=`echo $(printf "%0${nDigits}d" $1)`
else
	id=`echo $(printf "%0${nDigits}d" $SGE_TASK_ID)`
fi

method=astral
#echo working on $id
#mkdir -p $H/$id/$method
#cd $H/$id/$method
cd $H/$id

for i in $(find . -regex ".*\/.*\.tree[0-9]*")
#Removing weird errors
#rm -f original_g_trees.tree
#rm -f g_trees*.time
#cat g_trees* >> original_g_trees.tree
#done
#for i in $(find . -regex ".*\/original.*.tree")
do
	file=$(echo $i | sed -e "s/\.\/\(.*\)\.tree\(.*\)/\1.\2/g" -e "s/\.$//")
	leaves=($(cat $i |sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e 's/ /\'$'\n''/g' |sort|uniq|tail -n+2))
	species=($(echo ${leaves[*]}| tr ' ' '\n' |sed "s/\(.*\)\_.*\_.*$/\1/" | sort | uniq))
	ntaxa=${#species[@]}
	nloci=$(cat $i | wc -l)
	smap=""
	for sp in ${species[@]}
	do
        	individuals=($(echo ${leaves[*]} | tr ' ' '\n' | sed -n "/^${sp}_.*_.*/p"))
        	nind=${#individuals[*]}
        	smap=$(echo "$smap$sp $nind ${individuals[*]}\n")
	done
	echo -e $smap > ${file}.${method}.mapping

        if [[ ! -f ${file}.astral ]]
        then
		/usr/bin/time -a -o ${file}.${method}.time astral.sh -t 0 -s 2222 -a ${file}.${method}.mapping -i $i -o ${file}.astral 1> ${file}.${method}.out 2>${file}.${method}.err
        fi
done
