#!/bin/bash

dir="uly_retrieve"
# dir="tm_retrieve"

cd $dir
for atom in $(ls); do
    cd $atom
    for molecule in $(ls); do
        cd $molecule 
        for basis in $(ls); do
            cd $basis
            grep "TIMESTAMP" pyscf_output_mom.txt >> "temp_trip_stamp.txt"
            cd ..
        done
        cd ..
    done
    cd ..
done