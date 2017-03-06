#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

BIN=/home/dmallo/njstM/scripts/
DATA=/home/dmallo/njstM/data/broadsims/
nDigits=5

if [ "$SGE_TASK_ID" == "" ] || [ "$SGE_TASK_ID" == "undefined" ]
then
        id=`echo $(printf "%0${nDigits}d" $1)`
else
        id=`echo $(printf "%0${nDigits}d" $SGE_TASK_ID)`
fi

$BIN/get_time_folder.sh $DATA/$id
