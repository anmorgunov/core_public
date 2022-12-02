#!/bin/bash

tmpath="telemachus:/home/morgunov/core/calculations/mom/"
ulypath="ulysses:/work/morgunov/core/calculations/mom/"

uly_atoms="o"
# tm_atoms="f n"

mkdir -p "uly_retrieve"
mkdir -p "tm_retrieve"

# for uly_atom in $uly_atoms; do
#     scp -r "${ulypath}${uly_atom}/" "uly_retrieve/"
# done

# for tm_atom in $tm_atoms; do
#    scp -r "${tmpath}${tm_atom}/" "tm_retrieve/"
# done

rsync -avhu --progress uly_retrieve/ mom/
# rsync -avhu --progress tm_retrieve/ mom/
