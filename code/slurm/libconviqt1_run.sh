#!/bin/bash
#SBATCH --nodes=1
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

source $HOME/.bashrc
setup

srun -n 1 -c 64 $HOME/dependencies/bin/python3.7 ./../libconviqt1.py