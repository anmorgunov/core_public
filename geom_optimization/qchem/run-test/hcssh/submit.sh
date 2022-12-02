#!/bin/sh

#SBATCH -t 100:00:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=64000

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

echo `date` >> "time.txt"
qchem.latest qchem.in qchem.out
echo `date` >> "time.txt"
