from .config import SideChain
from .hierarchies import hierarchy_dict
import Bio.PDB as pdb
import numpy as np
from copy import deepcopy


class HydrogenBonding:
    _hydrogen_data = {"ALA": {},
                      "ARG": {"DONOR": {"NE": 1, "NH1": 2, "NH2": 2}},
                      "ASN": {"DONOR": {"ND2": 2},
                              "ACCEPTOR": {"OD1": 2}},
                      "ASP": {"ACCEPTOR": {"OD1": 2, "OD2": 2}},
                      "CYS": {},
                      "GLN": {"DONOR": {"NE2": 2},
                              "ACCEPTOR": {"OE1": 2}},
                      "GLU": {"ACCEPTOR": {"OE1": 2, "OE2": 2}},
                      "GLY": {},
                      "HIS": {"DONOR": {"ND1": 1, "NE2": 1},
                              "ACCEPTOR": {"ND1": 1, "NE2": 1}},
                      "ILE": {},
                      "LEU": {},
                      "LYS": {"DONOR": {"NZ": 3}},
                      "MET": {},
                      "PHE": {},
                      "PRO": {},
                      "SER": {"DONOR": {"OG": 1},
                              "ACCEPTOR": {"OG": 2}},
                      "THR": {"DONOR": {"OG1":1},
                              "ACCEPTOR": {"OG1": 2}},
                      "TRP": {"DONOR": {"NE1": 1}},
                      "TYR": {"DONOR": {"OH": 1},
                              "ACCEPTOR": {"OH": 1}},
                      "VAL": {}
                      }

    def __init__(self):
        pass

    def calc_existing_sidechain_hydrogen_bonds(self, chain, verbose=False):
        bonds = np.empty((0, 6))
        for i in range(len(chain) - 1):
            if len(self._hydrogen_data[chain[i].resname]) == 0:
                # cannot make hydrogen bond
                if verbose:
                    print(f"{i} ({chain[i].resname}) cannot make hydrogen bonds")
                continue

            i_at_hier = hierarchy_dict[chain[i].resname]
            i_at_hier = {k: [at for at in v if at[0] != 'H'] for k, v in i_at_hier.items()}
            i_roles = list(self._hydrogen_data[chain[i].resname].keys())
            if verbose:
                print(f"{i} ({chain[i].resname}) can be {i_roles}")

            for j in range(i + 1, len(chain)):
                if len(self._hydrogen_data[chain[j].resname]) == 0:
                    if verbose:
                        print(f"\t{j} ({chain[j].resname}) cannot make hydrogen bonds")
                    # cannot make hydrogen bond
                    continue
                j_roles = list(self._hydrogen_data[chain[j].resname].keys())
                if verbose:
                    print(f"\t{j} ({chain[j].resname}) can be {j_roles}")

                if len(set(i_roles + j_roles)) != 2:
                    # both are donors or acceptors
                    continue

                for i_role in i_roles:
                    for j_role in j_roles:
                        if i_role == j_role:
                            continue
                        for i_at_name in self._hydrogen_data[chain[i].resname][i_role]:
                            for j_at_name in self._hydrogen_data[chain[j].resname][j_role]:
                                res = calc_hydrogen_bond_likelihood(i_at_name,
                                                                    SideChain(i_at_hier,
                                                                              self._hydrogen_data,
                                                                              pdb_residue=chain[i]).get_atom_coords(),
                                                                    chain[i].resname,
                                                                    j_at_name,
                                                                    chain[j])
                                if verbose:
                                    print(f"\t\t{i_role}.{i_at_name} - {j_role}.{j_at_name}: {res}")
                                # print(res)
                                if res[0] < 0:
                                    if i_role == "DONOR":
                                        row = [i, chain[i].resname, i_at_name, j, chain[j].resname, j_at_name]
                                    else:
                                        row = [j, chain[j].resname, j_at_name, i, chain[i].resname, i_at_name]
                                    bonds = np.vstack((bonds, [row]))
        return bonds

    def get_available_hydrogen_bond_locations(self, existing_bonds):
        # existing_bonds = calc_existing_sidechain_hydrogen_bonds(chain)
        av_bonding = []
        for i in range(len(chain)):
            donor_pairs = np.where(existing_bonds[:, 0] == str(i))[0]
            acceptor_pairs = np.where(existing_bonds[:, 3] == str(i))[0]
            res_hydr_available = deepcopy(self._hydrogen_data[chain[i].resname])

            for k in donor_pairs:
                res_hydr_available["DONOR"][existing_bonds[k][2]] -= 1
                if res_hydr_available["DONOR"][existing_bonds[k][2]] < 0:
                    raise Exception("Wrong number of hydrogen bonds for single atom!")
            for k in acceptor_pairs:
                res_hydr_available["ACCEPTOR"][existing_bonds[k][5]] -= 1
                if res_hydr_available["ACCEPTOR"][existing_bonds[k][5]] < 0:
                    raise Exception("Wrong number of hydrogen bonds for single atom!")
            av_bonding += [res_hydr_available]
        return av_bonding

    def calc_hydrogen_bond_likelihood(self, i_at_name, tgt_at_coords, tgt_am_acid_type, j_at_name, j_res):
        # TODO: change to DSSP approach?

        d = np.linalg.norm(tgt_at_coords[i_at_name] - j_res[j_at_name].coord)
        if d > 3.5:
            return 0, "ERR_TOOFAR"
        i_roles = [k for k, v in self._hydrogen_data[tgt_am_acid_type].items() if i_at_name in v]
        j_roles = [k for k, v in self._hydrogen_data[j_res.resname].items() if j_at_name in v]
        if len(set(i_roles + j_roles)) != 2:
            return 0, "ERR_SINGLE_ROLE"
        if "DONOR" in i_roles and "ACCEPTOR" in j_roles:
            D = tgt_at_coords[i_at_name]
            A = j_res[j_at_name].coord
            X = tgt_at_coords[[k for k, v in hierarchy_dict[tgt_am_acid_type].items() if i_at_name in v][0]]
        elif "ACCEPTOR" in i_roles and "DONOR" in j_roles:
            A = tgt_at_coords[i_at_name]
            D = j_res[j_at_name].coord
            X = j_res[[k for k, v in hierarchy_dict[j_res.resname].items() if j_at_name in v][0]].coord
        else:
            return 0, "ERR_UNKNOWN_DONOR_ACCEPTOR_COMBINATION"

        ang = pdb.calc_angle(pdb.Vector(*X),
                             pdb.Vector(*D),
                             pdb.Vector(*A))
        if not np.pi / 2 <= ang <= np.pi:
            return 0, "ERR_WRONG_ANGLE"

        d0 = 3 / np.power(2, 1 / 6)

        #     return min((d0/d) ** 12 - (d0/d) ** 6,
        #                0)
        return min((d0 / d) ** 12 - (d0 / d) ** 6, 0), "OK"