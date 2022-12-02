#!/bin/bash

cur="$PWD/run/run-exp"
cd $cur
for mol in $(ls);
do
	cd $mol
	sbatch -J "geom-$mol" submit.sh
	cd ..
done
