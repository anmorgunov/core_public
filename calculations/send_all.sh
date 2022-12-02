#!/bin/bash

tmpath="telemachus:/home/morgunov/core/calculations"
ulypath="ulysses:/work/morgunov/core/calculations"
# paths="$tmpath $ulypath"
# paths="$tmpath"
paths="$ulypath"
files="run_mom.py"
folders="configs mols"

# for server_path in $paths; do
#     for file in $files; do
#         scp "$file" "$server_path/"
#     done

#     for folder in $folders; do
#         scp -r $folder $server_path
#     done
# done

# scp -r "submitscripts/telemachus" "$tmpath/submitscripts"
# scp "run_telemachus.sh" $tmpath
scp -r "submitscripts/ulysses" "$ulypath/submitscripts"
scp "run_ulysses.sh" $ulypath
