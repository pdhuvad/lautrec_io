#PBS -q  regular
#PBS -l mppwidth=256
#PBS -l walltime=02:00:00
#PBS -N D001-PBESol
#PBS -j eo

cd $PBS_O_WORKDIR

 echo "Chant and Be Happy" | mutt -s '2B2c-002-correct-lattice-has-started' 8595198518@tmomail.net
/usr/bin/python main.py
 echo "Chant and Be Happy" | mutt -s '2B2c-002-correct-lattice-has-ended' 8595198518@tmomail.net
