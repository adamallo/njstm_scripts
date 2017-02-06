#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH

module load python/2.7.8
H=/home/dmallo/njstM/data/broadsims/
BIN=/home/dmallo/njstM/scripts/

ASTRID="python /home/dmallo/astridm/ASTRID"
nDigits=5

if [ "$SGE_TASK_ID" == "" ] || [ "$SGE_TASK_ID" == "undefined" ]
then
	id=`echo $(printf "%0${nDigits}d" $1)`
else
	id=`echo $(printf "%0${nDigits}d" $SGE_TASK_ID)`
fi

#echo working on $id
#mkdir -p $H/$id/$method
#cd $H/$id/$method
cd $H/$id

for i in $(find . -regex ".*\/r.*\.tree[0-9]*")
do
	file=$(echo $i | sed -e "s/\.\/\(.*\)\.tree\(.*\)/\1.\2/g" -e "s/\.$//")
        cat $i |sed -e "s/:[^),]*//g" -e "s/)[0-9.]*//g" -e "s/[(,);]/ /g" -e 's/ /\'$'\n''/g' |sort|uniq|tail -n+2|sed "s/\(.*\)\_.*\_.*$/& \1/" > ${file}.astrid.mapping
        if [[ ! -f ${i}_basaltrifurcation ]]
        then
        	python $BIN/strictunroot.py -i $i -o ${i}_basaltrifurcation
        fi
        if [[ ! -f ${file}.astridmo ]]
        then
		method=astridmo
		/usr/bin/time -a -o ${file}.${method}.time $ASTRID -i ${i}_basaltrifurcation --map ${file}.astrid.mapping -c ${file}.${method}.cache -o ${file}.${method} 1> ${file}.${method}.out 2>${file}.${method}.err
        fi
	if [[ ! -f ${file}astridmu ]]
	then
		method=astridmu
		/usr/bin/time -a -o ${file}.${method}.time $ASTRID -i ${i}_basaltrifurcation --map ${file}.astrid.mapping -c ${file}.${method}.cache --bygene -o ${file}.${method} 1> ${file}.${method}.out 2>${file}.${method}.err
	fi
done
