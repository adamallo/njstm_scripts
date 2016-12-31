#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH

if [[ $# -ne 2 ]] || [[ ! -d $1 ]]
then
	echo "USAGE: $0 dir \"regexp\""
	echo "Regexp indicates the files to adapt after the leaf subsampling. For example: r*.tree"
        exit
fi

dir=$1
regexp=$2

for i in $dir/*
do
	if [[ -d $i ]]
	then
		for treefile in $i/${regexp}
		do
			if [[ ! -f ${treefile}.bkp ]]
			then
				cat $treefile | sed -e "s/^\[\&R\] //g" -e "s/'//g" > ${treefile}.temp
				mv $treefile ${treefile}.bkp
				mv ${treefile}.temp $treefile
			fi
		done
	fi
done
