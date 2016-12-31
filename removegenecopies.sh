#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH
#$ -V
module load python/2.7.8
python removegenecopies.py ../data/sim1/ 0.10 r0_10.tree -s 22
python removegenecopies.py ../data/sim1/ 0.25 r0_25.tree -s 22
python removegenecopies.py ../data/sim1/ 0.5 r0_50.tree -s 22
python removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 ../data/sim1/ 0.10 r1_10.tree -s 22
python removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 ../data/sim1/ 0.25 r1_25.tree -s 22
python removegenecopies.py -mk byindividual -ist 0.2 -itmin 0 -itmax 1 ../data/sim1/ 0.5 r1_50.tree -s 22
./simphygtrees_collapse.sh ../data/sim1/
