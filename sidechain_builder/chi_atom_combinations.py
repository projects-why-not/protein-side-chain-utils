
_atoms = {
    "ARG": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD"
      ],
      [
        "CB",
        "CG",
        "CD",
        "NE"
      ],
      [
        "CG",
        "CD",
        "NE",
        "CZ"
      ],
      [
        "CD",
        "NE",
        "CZ",
        "NH1"
      ]
    ],
    "ASN": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "OD1"
      ]
    ],
    "ASP": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "OD1"
      ]
    ],
    "CYS": [
      [
        "N",
        "CA",
        "CB",
        "SG"
      ]
    ],
    "GLN": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD"
      ],
      [
        "CB",
        "CG",
        "CD",
        "OE1"
      ]
    ],
    "GLU": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD"
      ],
      [
        "CB",
        "CG",
        "CD",
        "OE1"
      ]
    ],
    "HIS": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "ND1"
      ]
    ],
    "ILE": [
      [
        "N",
        "CA",
        "CB",
        "CG1"
      ],
      [
        "CA",
        "CB",
        "CG1",
        "CD1"
      ]
    ],
    "LEU": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD1"
      ]
    ],
    "LYS": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD"
      ],
      [
        "CB",
        "CG",
        "CD",
        "CE"
      ],
      [
        "CG",
        "CD",
        "CE",
        "NZ"
      ]
    ],
    "MET": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "SD"
      ],
      [
        "CB",
        "CG",
        "SD",
        "CE"
      ]
    ],
    "PHE": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD1"
      ]
    ],
    "PRO": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD"
      ]
    ],
    "SER": [
      [
        "N",
        "CA",
        "CB",
        "OG"
      ]
    ],
    "THR": [
      [
        "N",
        "CA",
        "CB",
        "OG1"
      ]
    ],
    "TRP": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD1"
      ]
    ],
    "TYR": [
      [
        "N",
        "CA",
        "CB",
        "CG"
      ],
      [
        "CA",
        "CB",
        "CG",
        "CD1"
      ]
    ],
    "VAL": [
      [
        "N",
        "CA",
        "CB",
        "CG1"
      ]
    ],
    "ALA": [],
    "GLY": []
}


def get_chi_atoms(resname, ang_num):
    return _atoms[resname][ang_num - 1]


def get_amino_acid_chi_num(resname):
    return len(_atoms[resname])
