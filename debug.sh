#PBS -q debug 
#PBS -l mppwidth=256
#PBS -l walltime=00:30:00
#PBS -N D001-PBESol
#PBS -j eo

cd $PBS_O_WORKDIR

/usr/bin/python main.py
