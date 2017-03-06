##I am not sure yet if I am going to use this data

qsub -q fast.q conduct_simulations.sh
qsub -q compute-1-x.q -pe threaded 8 scripts/removegenecopies.sh data/sim1/

./simphygtrees_collapse.sh ../data/sim1/

qsub -q compute-1-x.q scripts/post_removegenecopies.sh ".*\/r.*tree[0-9]*"

qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-49999 scripts/launch_njst.sh data/sim1/
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 50000-100000 scripts/launch_njst.sh data/sim1/

Rscript scripts/RF.r data/sim1/


##Broadsims
##Simulations from broadsims paper
##Scripts with H and ndigits were modified for this data (different directory and different length of the folder name)
qsub -q compute-1-x.q -pe threaded 8 scripts/removegenecopies.sh data/broadsims/
./simphygtrees_collapse.sh ../data/broadsims/
qsub -q compute-1-x.q scripts/post_removegenecopies.sh ".*\/r.*tree[0-9]*"

qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/launch_njst.sh data/broadsims/
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_astral.sh
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-10000 scripts/run_astridm.sh

Rscript scripts/RF.r data/broadsims/
