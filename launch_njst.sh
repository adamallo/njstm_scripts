#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

BIN=/home/dmallo/njstM/scripts/NJstM

if [ $# -lt 1 ]
then
        echo "USAGE: $0 [id] dir"
        exit
fi

if [ "$SGE_TASK_ID" == "" ] || [ "$SGE_TASK_ID" == "undefined" ]
then
        id=$1
	shift
else
        id=$SGE_TASK_ID
fi

#echo working on $id
dir=$1

module load R/3.2.2_1
id=$(printf "%05d" $id)

if [[ ! -d $dir/$id ]]
then
	echo "Error, the directory $dir/$id does not exist"
	echo "USAGE: $0 [id] dir"
	exit
else
	echo "Directory $dir/$id"
fi

cd $dir/$id
for i in $(find . -regex ".*\/r.*\.tree[0-9]*")
do
	file=$(echo $i | sed -e "s/\.\/\(.*\)\.tree\(.*\)/\1.\2/g" -e "s/\.$//")
	cat $i |sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e 's/ /\'$'\n''/g' |sort|uniq|tail -n+2|sed "s/\(.*\)\_.*\_.*$/& \1/" > ${file}.mapping
	if [[ ! -f ${file}.lnjst ]]
	then

		Rscript $BIN/njstm.r $i ${file}.mapping liu ${file}.lnjst
	fi
	if [[ ! -f ${file}.onjst ]]
	then

		Rscript $BIN/njstm.r $i ${file}.mapping original ${file}.onjst
	fi
	if [[ ! -f ${file}.unjst ]]
	then 
		Rscript $BIN/njstm.r $i ${file}.mapping unweighted ${file}.unjst
	fi
done
