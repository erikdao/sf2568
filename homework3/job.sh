#!/bin/bash -l

# time allocation
#SBATCH -A edu21.SF2568
# job name
#SBATCH -J myjob
# email notification
#SBATCH --mail-type=BEGIN,END
# 10 minutes wall-clock time will be given to this job
#SBATCH -t 00:10:00
# Number of nodes
#SBATCH --nodes=2
# set tasks per node to 24 in order to disable hyperthreading
#SBATCH --ntasks-per-node=24

module add i-compilers intelmpi

mpirun -np 48 ./hello

