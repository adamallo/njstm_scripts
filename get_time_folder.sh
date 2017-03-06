#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

rep=$1
cd $rep
echo "rep time mdata method" > time.stats
for j in *.time
do
	data=$(basename $j | sed "s/^\(.*\)\..*\.time/\1/g")
	method=$(basename $j | sed "s/^.*\.\(.*\)\.time/\1/g")
	time=$(tail -n 2 $j | awk 'BEGIN{FS=" "}{if (NR==1){sub("user","",$1); sub("user","",$2);print $1+$2}}')
	echo $rep $time $data $method >> time.stats
done
