#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

BIN=/home/dmallo/njstM/scripts/
H=/home/dmallo/njstM/data/broadsims/

echo "rep time mdata method" > time.stats
for i in $H/*
do
	if [[ -d $i ]]
	then
		rep=$(basename $i)
		echo Working in $i
		for j in $i/*.time
		do
			data=$(basename $j | sed "s/^\(.*\)\..*\.time/\1/g")
			method=$(basename $j | sed "s/^.*\.\(.*\)\.time/\1/g")
			time=$(tail -n 2 $j | awk 'BEGIN{FS=" "}{if (NR==1){sub("user","",$1); sub("user","",$2);print $1+$2}}')
			echo $rep $time $data $method >> time.stats
		done
	fi
done
