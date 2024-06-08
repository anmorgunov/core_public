import sys

sys.path.append(
    "../core_excitations"
)  # this allows us to import files from parent directories
import pprint

import tools

import PRIVATE

SUBMODULE = "basis_sets"  # how the current folder is named
BASE = PRIVATE.BASE_PATH + SUBMODULE


def parse_blocks(block1, block2):
    col1 = block1.split()
    cols = []
    for rline in block2.split("\n"):
        line = rline.split("\n")[0].split()
        nums = []
        for num in line:
            nums.append("{:e}".format(float(num)))
        cols.append(nums)
    final = ""
    for col1elt, othercols in zip(col1, cols):
        final += "           "
        final += "{:e}".format(float(col1elt))
        for otherelt in othercols:
            final += "           "
            final += otherelt
        final += "\n"
    # f = open(BASE+name, 'w')
    # f.write(final)
    return final


def parse_cfour(fname):

    with open(tools.join_path(BASE + [SUBMODULE, "cfour", fname]), "r") as f:
        cfour = f.read().split("\n\n")
    dicbloc = {
        "S": [cfour[1], cfour[2]],
        "P": [cfour[3], cfour[4]],
    }
    if fname != "D-H.cf":
        dicbloc["D"] = [cfour[5], cfour[6]]
    if fname != "T-H.cf" and "T" in fname.split("-")[0]:
        dicbloc["F"] = [cfour[7], cfour[8]]
    return dicbloc


def create_parsed_dic(element, dicbloc):
    outDic = {}
    for polarization, (block1, block2) in dicbloc.items():
        outDic[polarization] = parse_blocks(block1, block2)
    return outDic


def export_element(element, parsedic, zeta, isCore):
    zetaToPolToNum = {"D": {"S": 1, "P": 1}, "T": {"S": 2, "P": 2, "D": 1}}
    out = """"""
    for polarization, coeffs in parsedic.items():
        # print(coeffs)
        if isCore and polarization in zetaToPolToNum[zeta]:
            exclude = zetaToPolToNum[zeta][polarization]
            n_coeff = "\n".join(coeffs.split("\n")[: (-exclude - 1)])
            while True:
                array2D = [coeff.split() for coeff in n_coeff.split("\n")]
                if all(coeff[-1] == "0.000000e+00" for coeff in array2D):
                    no_zeros_array = [coeff[:-1] for coeff in array2D]
                    # print(no_zeros_array)
                    no_zeros = "\n".join(
                        ["           ".join(coeff) for coeff in no_zeros_array]
                    )
                    n_coeff = no_zeros
                else:
                    break
            # print(no_zeros_array)
            # print(no_zeros)
            out += f"{element} {polarization}\n{n_coeff}\n"
            for i in range(exclude):
                i_coeff = coeffs.split("\n")[-(i + 2)]
                no_zeros = "           ".join(
                    [coeff for coeff in i_coeff.split() if coeff != "0.000000e+00"]
                )
                # print(i_coeff)
                # print(no_zeros)
                out += f"{element} {polarization}\n{no_zeros}\n"
        else:
            out += f"{element} {polarization}\n{coeffs}\n"
    # print(out)
    return out


def compose_basis_set(zeta):
    atoms = "H C N O F"
    # atoms = 'N'
    basisSet = {}
    for atom in atoms.split():
        basisSet[atom] = {}
        dicbloc = parse_cfour(f"{zeta}-{atom}.cf")
        parsedic = create_parsed_dic(atom, dicbloc)
        basisSet[atom]["reg"] = export_element(atom, parsedic, zeta, isCore=False)
        if atom != "H":
            dicbloc = parse_cfour(f"C{zeta}-{atom}.cf")
            parsedic = create_parsed_dic(atom, dicbloc)
            basisSet[atom]["core"] = export_element(atom, parsedic, zeta, isCore=True)
    with open(
        tools.join_path(BASE + [SUBMODULE, "nwchem", f"cheng{zeta}.py"]), "w"
    ) as f:
        f.write("BASIS_SET = ")
        # f.write(json.dumps(basisSet))
        f.write(pprint.pformat(basisSet, width=1000))
    return basisSet


if __name__ == "__main__":
    compose_basis_set("D")
    compose_basis_set("T")

# Scheme
# For TZ: 2S, 2P, 1D
# For DZ: 1S, 1P

# dicbloc = parse_cfour('D-C.cf')
# o = create_parsed_dic('C', dicbloc)
# t = export_element('C', o)
# print(t)

# parse_cfour('T-C.cf')


# parse_s(s1, s2, 'cc-pVDZ_C-S.txt')
# parse_s(p1, p2, 'cc-pVDZ_C-P.txt')
# parse_s(d1, d2, 'cc-pVDZ_C-D.txt')
