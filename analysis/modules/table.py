import sys, constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools


class LaTeXTable:

    def __init__(self):
        pass

    def export_table(self, caption, label, positioning, headers, body, name):
        f = open(f"{name}.tex", 'w')
        f.write('\\begin{table}\n')
        f.write('  \caption{\\textbf{%s}}\n' % caption)
        f.write('  \label{tbl:%s}\n' % label)
        f.write('  \\begin{tabular}{%s}\n' % positioning)
        f.write('    \hline\n')
        printrow = lambda row: ['    ']+[f"{rowElt} & " if i+1 != len(row) else f"{rowElt} \\\ \n" for i, rowElt in enumerate(row)]
        if type(headers[0]) == list:
            for header in headers:
                f.write(''.join(printrow(header)))
            assert len(positioning.split()) == len(header), 'You must specify the alignment of each column'
        else:
            assert len(positioning.split()) == len(headers), 'You must specify the alignment of each column'
            f.write(''.join(printrow(headers)))
        f.write('    \hline\n')
        for row in body:
            f.write(''.join(printrow(row)))
        f.write('    \hline\n')
        f.write('  \end{tabular}\n')
        f.write('\end{table}\n')
        f.close()

if __name__ == "__main__":
    testHead = ['Molecule', 'Atom', 'UHF', 'UMP2']
    testBody = [['co', 'o', '12.3', '12.1'], ['co', 'c', '8.5', '8.4']]
    T = LaTeXTable()
    testCaption = 'This is an example of a table'
    testLabel = 'sampletable'
    testPositioning = 'c c c c'
    testPath = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'testTable'])
    T.export_table(testCaption, testLabel, testPositioning, testHead, testBody, testPath)