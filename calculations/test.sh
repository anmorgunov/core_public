#!/bin/bash

# old api
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/co.xyz" corebasis="cc-pCVDZ" regularbasis="cc-pVDZ" config=TEST
# new api 
python run_mom.py --element O --orbital 0 --geom "mols/exp/o/co.xyz" --corebasis "cc-pCVDZ" --regularbasis "cc-pVDZ" --config "test"

python run_mom.py --element C --orbital 1 --geom "mols/exp/c/c-o.xyz" --dospecialbasis --corebasis "pcX-1" --regularbasis "pc-1" --diis_start 5

python3 /home/am3939/core/calculations/run_mom.py atom=C orbital=1 geom=/home/am3939/core/calculations/mols/exp/c/c-o.xyz corebasis=pcX-1 regularbasis=pc-1 localize=false diis_start=9 doSpecialBasis=true > pyscf_output_mom.txt 
