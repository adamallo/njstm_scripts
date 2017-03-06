##ATTENTION: This script has been put together as a reference and using it directly would not necessarirly work properly due to directory changes etc.

qsub -q compute-1-x.q -pe threaded 8 scripts/removegenecopies.sh data/broadsims/

scripts/simphygtrees_collapse.sh data/broadsims/

qsub -q compute-1-x.q scripts/post_removegenecopies.sh ".*\/r.*tree[0-9]*"

qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/launch_njst.sh data/broadsims/
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_astral.sh
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_astridm.sh
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_time.sh
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_rf.sh

./scripts/concat_time.sh
./scripts/concat_rf.sh
./scripts/matrix.sh

echo "rep,nmiss,mdata";for j in *; do if [[ -d $j ]]; then id=$(basename $j); for i in $j/*.astridmu.cache; do nmiss=$(grep -c "\-99" $i); name=$(basename $i | sed "s/\.astridmu\.cache//g");echo $id,$nmiss,$name;done;fi;done > isc.csv

Rscript scripts/analyses.R
Rscript scripts/matrices.R
##OLD SIMULATIONS
#qsub -q fast.q conduct_simulations.sh
#qsub -q compute-1-x.q -pe threaded 8 scripts/removegenecopies.sh data/sim1/
#
#./simphygtrees_collapse.sh ../data/sim1/
#
#qsub -q compute-1-x.q scripts/post_removegenecopies.sh ".*\/r.*tree[0-9]*"
#
#qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-49999 scripts/launch_njst.sh data/sim1/
#qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 50000-100000 scripts/launch_njst.sh data/sim1/
#
#Rscript scripts/RF.r data/sim1/


