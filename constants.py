
ATOM_TO_MOLS = {
    'O': 'co co2'
}
from pprint import pprint 

FNAME_TO_MOLS = {
 'c-o': {'atom': 'C', 'latex': '\\textbf{C}O'},
 'c-o2': {'atom': 'C', 'latex': '\\textbf{C}O2'},
}

COLORS = {'red': ['#D7263D', '#A40E4C'], 
          'markers': ['#01BAEF',], # '#F18F01', '#4FB0C6'
          'line': ['#744FC6', '#D7BCE8']}

if __name__ == "__main__":
    d = {}
    for k, v in ATOM_TO_MOLS.items():
        d[k] = len(v.split())
    print(d)