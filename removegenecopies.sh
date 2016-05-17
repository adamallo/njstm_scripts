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
./simphygtrees_collapse.sh ../data/sim1/
