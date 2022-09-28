

# TODO: change order of hydrogens?


ala = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB1", "HB2", "HB3"]}

gly = {"N": ["CA"],
       "CA": ["HA2", "HA3"]}

val = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["CG1", "HB", "CG2"],
       "CG1": ["HG11", "HG12", "HG13"],
       "CG2": ["HG21", "HG22", "HG23"]}

leu = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["CD2", "HG", "CD1"],
       "CD1": ["HD11, HD12, HD13"],
       "CD2": ["HD21, HD22, HD23"]
       }

ile = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB", "CG1", "CG2"],
       "CG1": ["HG13", "CD1", "HG12"],
       "CG2": ["HG21", "HG22", "HG23"],
       "CD1": ["HD11", "HD12", "HD13"]}

met = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "SD", "HG2"],
       "SD": ["CE"],
       "CE": ["HE1", "HE2", "HE3"]}

phe = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["CD1", "CD2"],
       "CD1": ["CE1", "HD1"],
       "CD2": ["CE2", "HD2"],
       "CE1": ["CZ", "HE1"],
       "CE2": ["HE2"],
       "CZ": ["HZ"]}

trp = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["CD1", "CD2"],
       "CD1": ["NE1", "HD1"],
       "NE1": ["HE1"],
       "CD2": ["CE2", "CE3"],
       "CE2": ["CZ2"],
       "CZ2": ["CH2", "HZ2"],
       "CE3": ["CZ3", "HE3"],
       "CZ3": ["HZ3"],
       "CH2": ["HH2"]}

pro = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "CD", "HG2"],
       "CD": ["HD2", "HD3"]}

ser = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "OG", "HB3"],
       "OG": ["HG"]}

thr = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["OG1", "CG2", "HB"],
       "OG1": ["HG1"],
       "CG2": ["HG21", "HG22", "HG23"]}

cys = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "SG", "HB3"]}

tyr = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["CD1", "CD2"],
       "CD1": ["CE1", "HD1"],
       "CD2": ["CE2", "HD2"],
       "CE1": ["CZ", "HE1"],
       "CE2": ["HE2"],
       "CZ": ["OH"],
       "OH": ["HH"]}

asn = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["OD1", "ND2"],
       "ND2": ["HD21", "HD22"]}

gln = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "CD", "HG2"],
       "CD": ["OE1", "NE2"],
       "NE2": ["HE21", "HE22"]}

asp = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["OD1", "OD2"],
       "OD2": ["HD2"]}

glu = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "CD", "HG2"],
       "CD": ["OE1", "OE2"],
       "OD2": ["HD2"]}

lys = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "CD", "HG2"],
       "CD": ["HD2", "CE", "HD3"],
       "CE": ["HE3", "NZ", "HE2"],
       "NZ": ["HZ1", "HZ2", "HZ3"]}

arg = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["HG3", "CD", "HG2"],
       "CD": ["HD2", "NE", "HD3"],
       "NE": ["CZ", "HE"],
       "CZ": ["NH1", "NH2"],
       "NH1": ["HH11", "HH12"],
       "NH2": ["HH21", "HH22"]}

his = {"N": ["CA"],
       "CA": ["CB"],
       "CB": ["HB2", "CG", "HB3"],
       "CG": ["ND1", "CD2"],
       "ND1": ["CE1", "HD1"],
       "CE1": ["HE1"],
       "CD2": ["NE2", "HD2"],
       "NE2": ["HE2"]}

hierarchy_dict = {"ALA": ala,
                  "CYS": cys,
                  "ASP": asp,
                  "ASN": asn,
                  "GLU": glu,
                  "GLN": gln,
                  "GLY": gly,
                  "HIS": his,
                  "ILE": ile,
                  "LEU": leu,
                  "LYS": lys,
                  "MET": met,
                  "PHE": phe,
                  "PRO": pro,
                  "ARG": arg,
                  "SER": ser,
                  "THR": thr,
                  "VAL": val,
                  "TRP": trp,
                  "TYR": tyr
                  }
