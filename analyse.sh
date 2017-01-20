#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -v LD_LIBRARY_PATH
#$ -v PATH
module load R/3.2.2_1
module load gcc/5.2.0
Rscript /home/dmallo/njstM/scripts/RF.r $@
