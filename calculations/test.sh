#!/bin/bash

python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="cc-pCVDZ" regularbasis="cc-pVDZ" config=TEST
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="pcX-1" regularbasis="pc-1" doSpecialBasis=true
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="pcX-3" regularbasis="pc-3" doSpecialBasis=true config=TEST
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="pcX-4" regularbasis="pc-4" doSpecialBasis=true config=TEST

python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="ccX-DZ" regularbasis="cc-pVDZ" doSpecialBasis=true config=TEST
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="ccX-TZ" regularbasis="cc-pVTZ" doSpecialBasis=true config=TEST
python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="ccX-QZ" regularbasis="cc-pVQZ" doSpecialBasis=true config=TEST

python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/h2o.xyz" corebasis="ccX-TZ" regularbasis="cc-pCVTZ" doSpecialBasis=true config=TEST

python run_mom.py atom=O orbital=0 localize=False geom="mols/exp/o/co2.xyz" corebasis="pcX-1" regularbasis="pc-1" doSpecialBasis=true localize=true