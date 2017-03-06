#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH
DATA=/home/dmallo/njstM/data/broadsims/
nDigits=5
first=0
for i in $DATA/*
do
	if [[ -d $i ]]
	then
		if [[ $first -eq 0 ]]
		then
			first=1
			cat $i/time.stats > $DATA/time.stats
		else
			tail -n+2 $i/time.stats >> $DATA/time.stats
		fi
	fi
done
