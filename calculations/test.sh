#!/bin/bash

# old api
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/co.xyz" corebasis="cc-pCVDZ" regularbasis="cc-pVDZ" config=TEST
# new api 
python run_mom.py --element O --orbital 0 --geom "mols/exp/o/co.xyz" --corebasis "cc-pCVDZ" --regularbasis "cc-pVDZ" --config "test"