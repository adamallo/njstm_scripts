#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH
#$ -V
home="/home/dmallo/njstM/data"
simphy -rs 100000 -rl f:100 -rg 1 -sl f:20 -sb f:0.000001 -su f:sb -si u:1,5 -sp u:0.1E6,3E6 -cs 22 -v 1 -od 1 -op 1 -oc 1 -o ${home}/sim1 > ${home}/logs/simphy.out 2>&1
