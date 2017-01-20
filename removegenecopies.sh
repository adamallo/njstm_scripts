#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH
#$ -V
BIN=/home/dmallo/njstM/scripts/

if [[ $# -ne 1 ]]
then
	echo "Usage: $0 directory"
	exit
fi
cd $1
module load python/2.7.8
python $BIN/removegenecopies.py . 0.10 r0_10.tree -s 22 -ncores $NSLOTS
python $BIN/removegenecopies.py . 0.25 r0_25.tree -s 22 -ncores $NSLOTS
python $BIN/removegenecopies.py . 0.5 r0_50.tree -s 22 -ncores $NSLOTS
python $BIN/removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 . 0.10 r1_10.tree -s 22 -ncores $NSLOTS
python $BIN/removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 . 0.25 r1_25.tree -s 22 -ncores $NSLOTS
python $BIN/removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 . 0.5 r1_50.tree -s 22 -ncores $NSLOTS
#./simphygtrees_collapse.sh ../data/sim1/
