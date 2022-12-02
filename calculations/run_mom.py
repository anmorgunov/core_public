import numpy as np
import sys
import nwchem # see basis_sets module to understand where it comes from
import configs
from pyscf import gto, lo
from dmp2 import sgscf
from pyscf.tools import cubegen
from datetime import datetime


class CEBECalculator:

    def __init__(self, PARAMS):
        #heavyatom, nswap, geom, core_basis, full_basis, doLoc, testMode
        self.HEAVYATOM = PARAMS['atom']
        self.NSWAP = PARAMS['orbital']
        self.GEOM = PARAMS['geom']
        self.CORE_BASIS = PARAMS['corebasis']
        self.FULL_BASIS = PARAMS['regularbasis']
        self.DOLOC = PARAMS.get('localize', False)
        self.CONFIGKEY = PARAMS.get('config', self.FULL_BASIS)
        self.DIIS_START = PARAMS.get('diis_start', 0)
        self.DIIS_SPACE = PARAMS.get('diis_space', 50)
        self.MAX_SCF_CYCLES = PARAMS.get('scf_maxcycles', 1000)
        self.LEVEL_SHIFT = PARAMS.get('level_shift', 0)
        self.ATOMS_TO_LOCALIZE = PARAMS.get('atoms_to_localize', self.HEAVYATOM).split()
        

        self.externalBasisLibrary = {
            'cc-pVTZ': nwchem.deconccPVTZ.BASIS_SET,
            '4-31G': nwchem.g431.BASIS_SET,
            '6-311+G(3df)': nwchem.g6311plus3df.BASIS_SET,
        }
        self.keyToConfig = {
            'cc-pVQZ': configs.setups.Q_CONFIG,
            'cc-pV5Z': configs.setups.PENTUPLE_CONFIG,
            'test': configs.setups.TEST_CONFIG,
            'mp3': configs.setups.MP3_CONFIG
            }

        self.GEN_CONFIG = configs.setups.GEN_CONFIG

        self.config = self.keyToConfig.get(self.CONFIGKEY, self.GEN_CONFIG)
        # if 'cc-p' in self.FULL_BASIS:
        #     zeta = self.FULL_BASIS.split('V')[1][0]
        #     if zeta in self.keyToConfig:
        #         self.config = self.keyToConfig[zeta]
        
        # if self.CONFIGKEY == 'test': 
        #     self.config = self.keyToConfig[self.CONFIGKEY]

        

        self.doTime = configs.parameters.doTimestamps
        self.DO_T1_DIAGN = configs.parameters.DO_T1_DIAGN
        self.VERBOSITY = configs.parameters.VERBOSITY

        if self.doTime:
            self.timefile = open("timestamps.txt","w")
        
        self.doMP2 = self.config.get('doMP2', self.GEN_CONFIG['doMP2'])
        self.doMP3 = self.config.get('doMP3', self.GEN_CONFIG['doMP3'])
        self.doCCSD = self.config.get('doCCSD', self.GEN_CONFIG['doCCSD'])
        self.doTriples = self.config.get('doTriples', self.GEN_CONFIG['doTriples'])
        self.doSFX = self.config.get('doSFX', self.GEN_CONFIG['doSFX'])
        self.toFreeze = self.config.get('toFreeze', self.GEN_CONFIG['toFreeze'])
        self.doChengBasis = self.config.get('doChengBasis', self.GEN_CONFIG['doChengBasis'])
        self.doSpecialBasis = self.config.get('doSpecialBasis', self.GEN_CONFIG['doSpecialBasis'])
        self.basisLibKey = self.config.get('basisLibKey', self.GEN_CONFIG['basisLibKey'])
        self.toPrintDensity = self.config.get('toPrintDensity', self.GEN_CONFIG['toPrintDensity'])

    def get_time_now(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    #redacted

if __name__ == "__main__":
    '''
    Command line arguments are:
        geom: (string) Path to .xyz file
        core_basis (string): Basis for core atom
        full_basis (string): Basis for other atoms
    '''
    strToBool = {'false': False, 'true': True}
    parameters = sys.argv
    assert parameters[0].endswith('.py'), 'The first parameter should be path to the .py launch file'
    paramDic = {}
    case_sensitive = {'corebasis', 'regularbasis', 'atom', 'geom', 'atoms_to_localize'}
    for i, param in enumerate(parameters):
        if i == 0: continue
        key, val = param.split('=')
        if key.lower() in case_sensitive:
            paramDic[key] = val
        else:
            paramDic[key.lower()] = val.lower()
    assert 'geom' in paramDic, 'You must provide path to geometry file (use geom=/)'
    assert 'atom' in paramDic, 'You must provide symbol of an element to which core basis will be assigned to (use atom=)'
    assert 'orbital' in paramDic, 'You must provide index (0-based indexing) of an orbital to be excited (use orbital=)'
    assert 'regularbasis' in paramDic, 'You must specify the basis set (use regularbasis=)'
    assert 'corebasis' in paramDic, 'You must specify the basis set for atom that will be excited (use corebasis=)'
    if 'localize' in paramDic:
        paramDic['localize'] = strToBool[paramDic['localize']]

    toInt = 'orbital scf_maxcycles diis_start diis_space'.split()
    for key in toInt:
        if key in paramDic:
            paramDic[key] = int(paramDic[key])
    toInt = 'level_shift'.split()
    for key in toInt:
        if key in paramDic:
            paramDic[key] = float(paramDic[key])
    megalodon = CEBECalculator(paramDic)
    # print(paramDic)
    megalodon.run() 
