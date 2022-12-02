# import sys
# sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import tools 
import os 
import PRIVATE

CALC = PRIVATE.BASE_PATH + ['calculations']
from pprint import pprint

def parse_timestamps(path):
    times = {}
    with open(tools.join_path(path), 'r') as f:
        lines = f.readlines()
    for i, rline in enumerate(lines):
        line = rline.rstrip()
        if line.startswith('- '):
            method = line.split('- ')[1]
            times.setdefault(method, {})['start'] = lines[i+1]
            if len(lines) <= i+2: return None
            times.setdefault(method, {})['end'] = lines[i+2]
    return times

def parse_triples(path):
    times = {}
    with open(tools.join_path(path), 'r') as f:
        lines = f.readlines()
    times['CCSD(T)'] = lines[0].split('at ')[1]
    times['UCCSD(T)'] = lines[1].split('at ')[1]
    return times

def create_new_timestamps(path, times):
    if 'raw_timestamps.txt' not in (files:=os.listdir(tools.join_path(path, endWithSlash=True))):
        os.rename(tools.join_path(path + ['timestamps.txt']), tools.join_path(path + ['raw_timestamps.txt']))
    filestr = ''''''
    for method in 'HF MP2 MP3 CCSD CCSD(T) UHF UMP2 UMP3 UCCSD UCCSD(T)'.split():
        if path[-1] == 'Q' and method in {'MP3', 'UMP3'}: continue
        if path[-1] == '5' and method not in {'HF', 'UHF', 'MP2', 'UMP2'}: continue
        filestr += f"- {method}\n"
        filestr += times[method]['start']
        filestr += times[method]['end']
    with open(tools.join_path(path + ['timestamps.txt']), 'w') as f:
        f.write(filestr)


def combine_timestamp_files(path):
    old = parse_timestamps(path + ['timestamps.txt'])
    if old is None: return
    triples = parse_triples(path + ['temp_trip_stamp.txt'])
    tripleEnd = old['CCSD']['end']
    uTripleEnd = old['UCCSD']['end']
    old['CCSD']['end'] = triples['CCSD(T)']
    old.setdefault('CCSD(T)', {})['start'] = triples['CCSD(T)']
    old['CCSD(T)']['end'] = tripleEnd
    old['UCCSD']['end'] = triples['UCCSD(T)']
    old.setdefault('UCCSD(T)', {})['start'] = triples['UCCSD(T)']
    old['UCCSD(T)']['end'] = uTripleEnd
    create_new_timestamps(path, old)


def traverse(cluster):
    for atom in os.listdir(tools.join_path(CALC + [cluster])):
        for molecule in os.listdir(tools.join_path(CALC + [cluster, atom])):
            for basis in os.listdir(tools.join_path(CALC + [cluster, atom, molecule])):
                if basis == '5': continue
                path = CALC + [cluster, atom, molecule, basis]
                print(atom, molecule, basis)
                if 'raw_timestamps.txt' in (files:=os.listdir(tools.join_path(path, endWithSlash=True))): 
                    continue
                combine_timestamp_files(path)


if __name__ == "__main__":
    cluster = 'uly_retrieve'
    # cluster = 'tm_retrieve'

    traverse(cluster)