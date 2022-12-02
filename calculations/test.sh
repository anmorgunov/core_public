#!/bin/bash

mol="co2"
zeta="Q"

# python run_mom.py atom=Cl orbital=4 localize=True atoms_to_localize="C Cl" geom=mols/exp/cl/${mol}.xyz corebasis="cc-pCV${zeta}Z" regularbasis="cc-pV${zeta}Z" config=TEST # > mom_${mol}.txt
# python run_mom.py atom=O orbital=0 localize=False geom=mols/exp/o/${mol}.xyz corebasis="cc-pCV${zeta}Z" regularbasis="cc-pV${zeta}Z" #config=TEST # > mom_${mol}.txt
# python run_mom.py atom=F orbital=0 localize=False geom=mols/exp/f/${mol}.xyz corebasis="cc-pCV${zeta}Z" regularbasis="cc-pV${zeta}Z" config=TEST #> mom_${mol}.txt

basis="3-21G"
python run_mom.py atom=O orbital=0 localize=True geom=mols/exp/o/${mol}.xyz corebasis=${basis} regularbasis=${basis} config=TEST #> mom_${mol}.txt
