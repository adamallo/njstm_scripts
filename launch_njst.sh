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
id=$(printf "%06d" $id)

if [[ ! -d $dir/$id ]]
then
	echo "Error, the directory $dir/$id does not exist"
	echo "USAGE: $0 [id] dir"
	exit
else
	echo "Directory $dir/$id"
fi

cd $dir/$id
for i in *.tree
do
	file=$(basename $i .tree)
	cat $i |sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e 's/ /\'$'\n''/g' |sort|uniq|tail -n+2|sed "s/\(.*\)\_.*\_.*$/& \1/" > ${file}.mapping
	Rscript $BIN/njstm.r $i ${file}.mapping liu ${file}.lnjst
	Rscript $BIN/njstm.r $i ${file}.mapping original ${file}.onjst
	Rscript $BIN/njstm.r $i ${file}.mapping unweighted ${file}.unjst
done
