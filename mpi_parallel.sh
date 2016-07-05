#!/bin/bash --login
#PBS -q short
#PBS -l nodes=4:ppn=16
#PBS -l walltime=4:00:00
#PBS -N 11-000
#PBS -j oe
#PBS -m n



cd ${PBS_O_WORKDIR}
#perl xifan_stress_full> log_run
#python pratik_psinv_gen_lautrec.py >log_run
#mpirun -np 20 --map-by socket ~/local/vasp.5.3/vasp

python main.py>log_run

#mkdir run1
#cp log* out* ./run1
#
#/usr/bin/python main.py
#
#
#mkdir run2
#cp log* out* ./run2
#
#/usr/bin/python main.py
#
#mkdir run3
#cp log* out* ./run3
#
#/usr/bin/python main.py
#mkdir run4
#cp log* out* ./run4
#
#/usr/bin/python main.py
#
