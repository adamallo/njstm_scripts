qsub -q fast.q conduct_simulations.sh

python removegenecopies.py ../data/sim1/ 0.10 r0_10.tree -s 22

./simphygtrees_collapse.sh ../data/sim1/

qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 1-49999 scripts/launch_njst.sh data/sim1/
qsub -q compute-1-x.q,compute-0-x.q -j y -e ./e_logs/ -o ./o_logs/ -t 50000-100000 scripts/launch_njst.sh data/sim1/

