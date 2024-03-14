from pprint import pprint

ATOM_TO_MOLS = {
    "C": "c-cl2h2 c-cl3h c-cl4 c-clh3 c-f2h2 c-f3h c-f4 c-fh3 c-h3cch c-h3cn ch3-c-n c-h3nc c-h3oh c-h4 c-o2 h2c-c-ch2 h-c-o2ch3 c2-h2 c2-h4 c2-h6 c2h6-c-o c-h2o ch3-c-ho ch3-c-o2h c-o c-h3och3",  # ch3n-c n-c-cn h-c-n
    "N": "c2h5nh2 ch3nh2 ch3nhch3 ch3x3n hcn hconh2 hno2 hno3 nfh2 nh3 nonh2 pyrrole n2h4 pyridine"
    + " c2h5c-n ch2chc-n ch3sc-n h-n-co h2n-c2h4oh h2n-cho hco-n-hch3 h2n-conh2 n-cch2cn"
    + " ofpyridine pfpyridine me2-n-cho ipr-nh2 h2n-c2h4nh2 pr-nh2 pohpyridine mnh2pyridi-n-e m-n-h2pyridine o-n-h2pyridine pnh2pyridi-n-e p-n-h2pyridine 2pyrido-n-e 2pyrido-n-earo onh2pyridi-n-e",  # + ' ch3-n-o2', # c2n2 #
    "O": "h2o c2h5oh c4h4o ch3cooh ch3coch3 ch3och3 ch3oh co hcho hcooch3 hcooh co2"
    + " cf3c-o-oh cf3co-o-h ch2chcho ch3c-o-oh ch3n-o2 h2nc-o-nh2 h2nch-o hc-o-och3 hnc-o ipr-oh pr-oh ch3c-o-och3 ch3co-o-ch3",  #
    "F": "c2h3f c2h5f c3h7f ch3cof ch3nhcof ch3nhf ch3of foh hccf hoch2f nh2ch2f nh2f ch3f chf3 ch2f2 cf4 f2 c5h5f"
    + " hf pf3 cf3chch2 cf3cch bf3 cf3ocf3",  # c2f4
}

FNAME_TO_MOLS = {
    "c-cl2h2": {
        "atom": "C",
        "latex": "\\textbf{C}H2Cl2",
        "formula": "CH2Cl2",
    },
    "c-cl3h": {
        "atom": "C",
        "latex": "\\textbf{C}HCl3",
        "formula": "CHCl3",
    },
    "c-cl4": {
        "atom": "C",
        "latex": "\\textbf{C}Cl4",
        "formula": "CCl4",
    },
    "c-clh3": {
        "atom": "C",
        "latex": "\\textbf{C}H3Cl",
        "formula": "CH3Cl",
    },
    "c-f2h2": {
        "atom": "C",
        "latex": "\\textbf{C}H2F2",
        "formula": "CH2F2",
    },
    "c-f3h": {
        "atom": "C",
        "latex": "\\textbf{C}HF3",
        "formula": "CHF3",
    },
    "c-f4": {
        "atom": "C",
        "latex": "\\textbf{C}F4",
        "formula": "CF4",
    },
    "c-fh3": {
        "atom": "C",
        "latex": "\\textbf{C}H3F",
        "formula": "CH3F",
    },
    "c-h2o": {
        "atom": "C",
        "latex": "H\\textbf{C}HO",
        "formula": "HCHO",
    },
    "c-h3cch": {
        "atom": "C",
        "latex": "\\textbf{C}H3CCH",
        "formula": "CH3CCH",
    },
    "c-h3cn": {
        "atom": "C",
        "latex": "\\textbf{C}H3CN",
        "formula": "CH3CN",
    },
    "c-h3nc": {
        "atom": "C",
        "latex": "\\textbf{C}H3NC",
        "formula": "CH3NC",
    },
    "c-h3oh": {
        "atom": "C",
        "latex": "\\textbf{C}H3OH",
        "formula": "CH3OH",
    },
    "c-h3och3": {
        "atom": "C",
        "latex": "\\textbf{C}H3OCH3",
        "formula": "CH3OCH3",
    },
    "c-h4": {
        "atom": "C",
        "latex": "\\textbf{C}H4",
        "formula": "CH4",
    },
    "c-o": {
        "atom": "C",
        "latex": "\\textbf{C}O",
        "formula": "CO",
    },
    "c-o2": {
        "atom": "C",
        "latex": "\\textbf{C}O2",
        "formula": "CO2",
    },
    "c2-h2": {
        "atom": "C",
        "latex": "\\textbf{C}2H2",
        "formula": "C2H2",
    },
    "c2-h4": {
        "atom": "C",
        "latex": "\\textbf{C}2H4",
        "formula": "C2H4",
    },
    "c2-h6": {
        "atom": "C",
        "latex": "\\textbf{C}2H6",
        "formula": "C2H6",
    },
    "c2h3f": {
        "atom": "F",
        "latex": "C2H3\\textbf{F}",
        "formula": "C2H3F",
    },
    "c2h5f": {
        "atom": "F",
        "latex": "C2H5\\textbf{F}",
        "formula": "C2H5F",
    },
    "c2h5nh2": {
        "atom": "N",
        "latex": "C2H5\\textbf{N}H2",
        "formula": "C2H5NH2",
    },
    "c2h5oh": {
        "atom": "O",
        "latex": "C2H5\\textbf{O}H",
        "formula": "C2H5OH",
    },
    "c2h6-c-o": {
        "atom": "C",
        "latex": "(CH3)2\\textbf{C}O",
        "formula": "(CH3)2CO",
    },
    "c3h7f": {
        "atom": "F",
        "latex": "C3H7\\textbf{F}",
        "formula": "C3H7F",
    },
    "c4h4o": {
        "atom": "O",
        "latex": "C4H4\\textbf{O}",
        "formula": "C4H4O",
    },
    "c5h5f": {
        "atom": "F",
        "latex": "C5H5\\textbf{F}",
        "formula": "C5H5F",
    },
    "cf4": {
        "atom": "F",
        "latex": "C\\textbf{F}4",
        "formula": "CF4",
    },
    "ch2f2": {
        "atom": "F",
        "latex": "CH2\\textbf{F}2",
        "formula": "CH2F2",
    },
    "ch3-c-ho": {
        "atom": "C",
        "latex": "CH3\\textbf{C}HO",
        "formula": "CH3CHO",
    },
    "ch3-c-n": {
        "atom": "C",
        "latex": "CH3\\textbf{C}N",
        "formula": "CH3CN",
    },
    "ch3-c-o2h": {
        "atom": "C",
        "latex": "CH3\\textbf{C}O2H",
        "formula": "CH3CO2H",
    },
    "ch3coch3": {
        "atom": "O",
        "latex": "(CH3)2C\\textbf{O}",
        "formula": "(CH3)2CO",
    },
    "ch3cof": {
        "atom": "F",
        "latex": "CH3CO\\textbf{F}",
        "formula": "CH3COF",
    },
    "ch3cooh": {
        "atom": "O",
        "latex": "CH3CO\\textbf{O}H",
        "formula": "CH3COOH",
    },
    "ch3f": {
        "atom": "F",
        "latex": "CH3\\textbf{F}",
        "formula": "CH3F",
    },
    "ch3nh2": {
        "atom": "N",
        "latex": "CH3\\textbf{N}H2",
        "formula": "CH3NH2",
    },
    "ch3nhch3": {
        "atom": "N",
        "latex": "CH3\\textbf{N}HCH3",
        "formula": "CH3NHCH3",
    },
    "ch3nhcof": {
        "atom": "F",
        "latex": "CH3NHCO\\textbf{F}",
        "formula": "CH3NHCOF",
    },
    "ch3nhf": {
        "atom": "F",
        "latex": "CH3NH\\textbf{F}",
        "formula": "CH3NHF",
    },
    "ch3och3": {
        "atom": "O",
        "latex": "CH3\\textbf{O}CH3",
        "formula": "CH3OCH3",
    },
    "ch3of": {
        "atom": "F",
        "latex": "CH3O\\textbf{F}",
        "formula": "CH3OF",
    },
    "ch3oh": {
        "atom": "O",
        "latex": "CH3\\textbf{O}H",
        "formula": "CH3OH",
    },
    "ch3x3n": {
        "atom": "N",
        "latex": "(CH3)3\\textbf{N}",
        "formula": "(CH3)3N",
    },
    "chf3": {
        "atom": "F",
        "latex": "CH\\textbf{F}3",
        "formula": "CHF3",
    },
    "co": {
        "atom": "O",
        "latex": "C\\textbf{O}",
        "formula": "CO",
    },
    "co2": {
        "atom": "O",
        "latex": "C\\textbf{O}2",
        "formula": "CO2",
    },
    "f2": {
        "atom": "F",
        "latex": "\\textbf{F}2",
        "formula": "F2",
    },
    "hf": {
        "atom": "F",
        "latex": "H\\textbf{F}",
        "formula": "HF",
    },
    "foh": {
        "atom": "F",
        "latex": "HO\\textbf{F}",
        "formula": "HOF",
    },
    "h-c-o2ch3": {
        "atom": "C",
        "latex": "H\\textbf{C}O2CH3",
        "formula": "HCO2CH3",
    },
    "h2c-c-ch2": {
        "atom": "C",
        "latex": "H2C\\textbf{C}CH2",
        "formula": "H2CCCH2",
    },
    "h2o": {
        "atom": "O",
        "latex": "H2\\textbf{O}",
        "formula": "H2O",
    },
    "hccf": {
        "atom": "F",
        "latex": "HCC\\textbf{F}",
        "formula": "HCCF",
    },
    "hcho": {
        "atom": "O",
        "latex": "HCH\\textbf{O}",
        "formula": "HCHO",
    },
    "hcn": {
        "atom": "N",
        "latex": "HC\\textbf{N}",
        "formula": "HCN",
    },
    "hconh2": {
        "atom": "N",
        "latex": "HCO\\textbf{N}H2",
        "formula": "HCONH2",
    },
    "hcooch3": {
        "atom": "O",
        "latex": "HCO\\textbf{O}CH3",
        "formula": "HCOOCH3",
    },
    "hcooh": {
        "atom": "O",
        "latex": "HCO\\textbf{O}H",
        "formula": "HCOOH",
    },
    "hno2": {
        "atom": "N",
        "latex": "H\\textbf{N}O2",
        "formula": "HNO2",
    },
    "hno3": {
        "atom": "N",
        "latex": "H\\textbf{N}O3",
        "formula": "HNO3",
    },
    "hoch2f": {
        "atom": "F",
        "latex": "HOCH2\\textbf{F}",
        "formula": "HOCH2F",
    },
    "n2h4": {
        "atom": "N",
        "latex": "\\textbf{N}2H4",
        "formula": "N2H4",
    },
    "nfh2": {
        "atom": "N",
        "latex": "\\textbf{N}H2F",
        "formula": "NH2F",
    },
    "nh2ch2f": {
        "atom": "F",
        "latex": "NH2CH2\\textbf{F}",
        "formula": "NH2CH2F",
    },
    "nh2f": {
        "atom": "F",
        "latex": "NH2\\textbf{F}",
        "formula": "NH2F",
    },
    "nh3": {
        "atom": "N",
        "latex": "\\textbf{N}H3",
        "formula": "NH3",
    },
    "nonh2": {
        "atom": "N",
        "latex": "NO\\textbf{N}H2",
        "formula": "NONH2",
    },
    "pyridine": {
        "atom": "N",
        "latex": "C5H5\\textbf{N}",
        "formula": "C5H5N",
    },
    "pyrrole": {
        "atom": "N",
        "latex": "C4H5\\textbf{N}",
        "formula": "C4H5N",
    },
    "bf3": {
        "atom": "F",
        "latex": "B\\textbf{F}3",
        "formula": "BF3",
    },
    "c2h5c-n": {
        "atom": "N",
        "latex": "C2H5C\\textbf{N}",
        "formula": "C2H5CN",
    },
    "cf3c-o-oh": {
        "atom": "O",
        "latex": "CF3C\\textbf{O}OH",
        "formula": "CF3COOH",
    },
    "cf3cch": {
        "atom": "F",
        "latex": "C\\textbf{F}3CCH",
        "formula": "CF3CCH",
    },
    "cf3chch2": {
        "atom": "F",
        "latex": "C\\textbf{F}3CHCH2",
        "formula": "CF3CHCH2",
    },
    "cf3co-o-h": {
        "atom": "O",
        "latex": "CF3CO\\textbf{O}H",
        "formula": "CF3COOH",
    },
    "ch2chc-n": {
        "atom": "N",
        "latex": "CH2CHC\\textbf{N}",
        "formula": "CH2CHCN",
    },
    "ch2chcho": {
        "atom": "O",
        "latex": "CH2CHCH\\textbf{O}",
        "formula": "CH2CHCHO",
    },
    "ch3-n-o2": {
        "atom": "N",
        "latex": "CH3\\textbf{N}O2",
        "formula": "CH3NO2",
    },
    "ch3c-o-oh": {
        "atom": "O",
        "latex": "CH3C\\textbf{O}OH",
        "formula": "CH3COOH",
    },
    "ch3n-o2": {
        "atom": "O",
        "latex": "CH3N\\textbf{O}2",
        "formula": "CH3NO2",
    },
    "ch3sc-n": {
        "atom": "N",
        "latex": "CH3SC\\textbf{N}",
        "formula": "CH3SCN",
    },
    "h-n-co": {
        "atom": "N",
        "latex": "H\\textbf{N}CO",
        "formula": "HNCO",
    },
    "h2n-c2h4oh": {
        "atom": "N",
        "latex": "H2\\textbf{N}C2H4OH",
        "formula": "H2NC2H4OH",
    },
    "h2n-cho": {
        "atom": "N",
        "latex": "H2\\textbf{N}CHO",
        "formula": "H2NCHO",
    },
    "h2n-conh2": {
        "atom": "N",
        "latex": "H2\\textbf{N}CONH2",
        "formula": "H2NCONH2",
    },
    "h2nc-o-nh2": {
        "atom": "O",
        "latex": "H2NC\\textbf{O}NH2",
        "formula": "H2NCONH2",
    },
    "h2nch-o": {
        "atom": "O",
        "latex": "H2NCH\\textbf{O}",
        "formula": "H2NCHO",
    },
    "hc-o-och3": {
        "atom": "O",
        "latex": "HC\\textbf{O}OCH3",
        "formula": "HCOOCH3",
    },
    "hco-n-hch3": {
        "atom": "N",
        "latex": "HCO\\textbf{N}HCH3",
        "formula": "HCONHCH3",
    },
    "hf": {
        "atom": "F",
        "latex": "H\\textbf{F}",
        "formula": "HF",
    },
    "hnc-o": {
        "atom": "O",
        "latex": "HNC\\textbf{O}",
        "formula": "HNCO",
    },
    "ipr-oh": {
        "atom": "O",
        "latex": "i-Pr\\textbf{O}H",
        "formula": "i-PrOH",
    },
    "n-cch2cn": {
        "atom": "N",
        "latex": "\\textbf{N}CCH2CN",
        "formula": "NCCH2CN",
    },
    "pf3": {
        "atom": "F",
        "latex": "P\\textbf{F}3",
        "formula": "PF3",
    },
    "pr-oh": {
        "atom": "O",
        "latex": "Pr\\textbf{O}H",
        "formula": "PrOH",
    },
    "2pyrido-n-e": {
        "atom": "N",
        "latex": "C5H5\\textbf{N}O",
        "formula": "C5H5NO",
    },
    "2pyrido-n-earo": {
        "atom": "N",
        "latex": "aro-C5H5\\textbf{N}O",
        "formula": "aro-C5H5NO",
    },
    "cf3ocf3": {
        "atom": "F",
        "latex": "C\\textbf{F}3OCF3",
        "formula": "CF3OCF3",
    },
    "ch3c-o-och3": {
        "atom": "O",
        "latex": "CH3C\\textbf{O}OCH3",
        "formula": "CH3COOCH3",
    },
    "ch3co-o-ch3": {
        "atom": "O",
        "latex": "CH3CO\\textbf{O}CH3",
        "formula": "CH3COOCH3",
    },
    "h2n-c2h4nh2": {
        "atom": "N",
        "latex": "H2\\textbf{N}C2H4NH2",
        "formula": "H2NC2H4NH2",
    },
    "ipr-nh2": {
        "atom": "N",
        "latex": "i-Pr\\textbf{N}H2",
        "formula": "i-PrNH2",
    },
    "m-n-h2pyridine": {
        "atom": "N",
        "latex": "m-\\textbf{N}H2-C5H4N",
        "formula": "m-NH2-C5H4N",
    },
    "me2-n-cho": {
        "atom": "N",
        "latex": "(CH3)2\\textbf{N}CHO",
        "formula": "(CH3)2NCHO",
    },
    "mnh2pyridi-n-e": {
        "atom": "N",
        "latex": "m-NH2-C5H4\\textbf{N}",
        "formula": "m-NH2-C5H4N",
    },
    "o-n-h2pyridine": {
        "atom": "N",
        "latex": "o-\\textbf{N}H2-C5H4N",
        "formula": "o-NH2-C5H4N",
    },
    "ofpyridine": {
        "atom": "N",
        "latex": "o-F-C5H4\\textbf{N}",
        "formula": "o-F-C5H4N",
    },
    "onh2pyridi-n-e": {
        "atom": "N",
        "latex": "o-NH2-C5H4\\textbf{N}",
        "formula": "o-NH2-C5H4N",
    },
    "p-n-h2pyridine": {
        "atom": "N",
        "latex": "p-\\textbf{N}H2-C5H4N",
        "formula": "p-NH2-C5H4N",
    },
    "pfpyridine": {
        "atom": "N",
        "latex": "p-F-C5H4\\textbf{N}",
        "formula": "p-F-C5H4N",
    },
    "pnh2pyridi-n-e": {
        "atom": "N",
        "latex": "p-NH2-C5H4\\textbf{N}",
        "formula": "p-NH2-C5H4N",
    },
    "pohpyridine": {
        "atom": "N",
        "latex": "p-OH-C5H4\\textbf{N}",
        "formula": "p-OH-C5H4N",
    },
    "pr-nh2": {"atom": "N", "latex": "Pr-\\textbf{N}H2", "formula": "Pr-NH2"},
}

COLORS = {
    "red": ["#D7263D", "#A40E4C"],
    "markers": [
        "#01BAEF",
    ],  # '#F18F01', '#4FB0C6'
    "line": ["#744FC6", "#D7BCE8"],
}

mp2_opts = "MP2(T), MP2(Q), MP2(5), MP2(T Q), MP2(D T Q), MP2(T Q 5), MP2(D T Q 5), MP2(D T Q 5)"
mp2corr_opts = "MP2(T)+DifD, MP2(T)+DifD(T), MP2(Q)+DifD, MP2(Q)+DifD(T), MP2(5)+DifD, MP2(5)+DifD(T), MP2(T Q)+DifD, MP2(T Q)+DifD(T), MP2(D T Q)+DifD, MP2(D T Q)+DifD(T), MP2(D T)+DifD, MP2(D T)+DifD(T), MP2(T Q 5)+DifD, MP2(T Q 5)+DifD(T), MP2(D T Q 5)+DifD, MP2(D T Q 5)+DifD(T)"
EXTRAPOLATION_OPTIONS = {
    "HF": "T-Q-HF D-T-Q-HF T-Q-5-HF D-T-Q-5-HF".split(),
    "CCSD": "D-T-CCSD D-T-CCSD(T) T-Q-CCSD T-Q-CCSD(T) D-T-Q-CCSD D-T-Q-CCSD(T)".split(),
    "MP2": mp2_opts.split(", "),
    "MP2+Corr": mp2corr_opts.split(", "),
    # + "".split(
    #     ", "
    # )
    # + "".split(
    #     ", "
    # )
    # + "".split(
    #     ", "
    # )
    # + "".split(", "),
    # "ccsdTriplesStudy": "T-Q-CCSD T-Q-CCSD(T)".split(),
    # "MP2triplesStudy": "MP2(Q)+DifD, MP2(Q)+DifD(T), MP2(5)+DifD, MP2(5)+DifD(T), MP2(T Q)+DifD, MP2(T Q)+DifD(T)".split(
    #     ", "
    # ),
    # "bestOptions": "T-Q-CCSD, T-Q-CCSD(T), MP2(Q)+DifD, MP2(T Q)+DifD, MP2(5)+DifD".split(
    #     ", "
    # ),
    # "MP2difStudy": "MP2(T)+DifD, MP2(T), MP2(Q)+DifD, MP2(Q), MP2(5)+DifD, MP2(5), MP2(T Q)+DifD, MP2(T Q), MP2(T Q 5)+DifD, MP2(T Q 5)".split(
    #     ", "
    # ),
}
# pprint(EXTRAPOLATION_OPTIONS)
# self.all = (
#     self.HFoptions + self.CCSDoptions + self.MP2options
# )  # + self.MP2triplesStudy + self.bestOptions

if __name__ == "__main__":
    d = {}
    for k, v in ATOM_TO_MOLS.items():
        d[k] = len(v.split())
    print(d)
