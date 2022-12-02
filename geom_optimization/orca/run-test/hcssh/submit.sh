#!/bin/sh

#SBATCH -t 100:00:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=64000

echo `date` >> "time.txt"
/home/morgunov/orca/orca geom.inp > geom.out
echo `date` >> "time.txt"
