from .advanced_geometry import get_axis_rotation_matrix, project_point_on_line
from .volume import VolumeProvider
from .hierarchies import hierarchy_dict
from .chi_atom_combinations import get_chi_atoms, get_amino_acid_chi_num
from .potential.rotameric import get_normalized_rot_potential
import numpy as np
from Bio.PDB import calc_angle, calc_dihedral, Vector


class AtomConfig:
    # _children = None
    # _parent = None
    _fix_coord = None

    def __init__(self, name, dihedral=None, planar=None, bond_len=None, fix_coord=None):
        self._name = name
        if fix_coord is not None and name in ["N", "CA", "CB"]:
            self._fix_coord = fix_coord
        else:
            self._dihedral = dihedral
            self._planar = planar
            self._bond_len = bond_len
            self._fix_coord = None

    def set_geometry(self, geom_triple=None, fix_coord=None):
        if geom_triple is None and fix_coord is None:
            raise Exception("No values provided to atom.set_geometry!")
        if fix_coord is not None:
            self._fix_coord = fix_coord
        elif geom_triple is not None:
            if len(geom_triple) != 3:
                raise Exception("Wrong geom_triple shape!")
            self._dihedral, self._planar, self._bond_len = geom_triple

    def set_dihedral(self, val):
        self._dihedral = val

    def get_dihedral(self):
        return self._dihedral

    # def set_planar(self, val):
    #     self._planar = val
    #
    # def set_bond_len(self, val):
    #     self._bond_len = val

    def get_coord(self, prev_atoms):
        # print(f"get coord: {self._name}, prev = {prev_atoms}")
        if self._fix_coord is not None:
            return self._fix_coord

        if prev_atoms.shape != (3,3):
            raise Exception("Wrong prev_atoms shape!")
        if self._dihedral is None or self._bond_len is None or self._planar is None:
            raise Exception("Atom geometrical properties not initialized!")
        return self.place_next_dihedral_atom(prev_atoms)

    def place_next_dihedral_atom(self, prev_atoms):
        A1 = prev_atoms[-2] + get_axis_rotation_matrix(prev_atoms[-1] - prev_atoms[-2],
                                                       self._dihedral).dot(prev_atoms[-3] - prev_atoms[-2])
        v0 = (prev_atoms[-1] - prev_atoms[-2]) / np.linalg.norm(prev_atoms[-1] - prev_atoms[-2])
        E = project_point_on_line(A1,
                                  (prev_atoms[-1],
                                   prev_atoms[-1] - prev_atoms[-2]))
        v1 = (A1 - E) / np.linalg.norm(A1 - E)
        return prev_atoms[-1] - v0 * self._bond_len * np.cos(self._planar) + v1 * self._bond_len * np.sin(self._planar)


class SideChain:
    def __init__(self,
                 resname,
                 hydrogen_bond_info,
                 pdb_residue=None,
                 drop_hydrogens=True):
        """
        :param atom_hierarchy_dict: dict; corresponding dict from ./hierarchies module.
        :param hydrogen_bond_info: dict: keys: "DONOR", "ACCEPTOR"; values: lists of corresponding atoms. If type
                is "NONE", this parameter should be equal to {}.
        :param pdb_residue: Bio.PDB.Residue; residue with original geometry.
        """
        self.resname = resname
        self._hierarchy = hierarchy_dict[resname]
        if drop_hydrogens:
            self._hierarchy = {k:[at for at in v if at[0] != "H"] for k,v in self._hierarchy.items()}

        self._hydrogen_bonding = hydrogen_bond_info
        atom_names = list(self._hierarchy.keys())
        atom_names += [item for subl in self._hierarchy.values() for item in subl]
        atom_names = list(set(atom_names))
        self._atoms = {}
        for at_name in atom_names:
            self._atoms[at_name] = AtomConfig(at_name)
        if pdb_residue is not None:
            self._init_w_pdb_residue(pdb_residue)

    def _init_w_pdb_residue(self, residue):
        self._init_branch_w_pdb_residue(["N"], residue)

    def _init_branch_w_pdb_residue(self, atom_chain, residue):
        # setting last (target) atom in current chain
        if len(atom_chain) < 4:
            # if last atom is N, CA, CB, init by original coord
            self._atoms[atom_chain[-1]].set_geometry(fix_coord=residue[atom_chain[-1]].coord)
        else:
            if residue.has_id(atom_chain[-1]):
                dih = calc_dihedral(residue[atom_chain[-4]].get_vector(),
                                    residue[atom_chain[-3]].get_vector(),
                                    residue[atom_chain[-2]].get_vector(),
                                    residue[atom_chain[-1]].get_vector()
                                    )
                planar = calc_angle(residue[atom_chain[-3]].get_vector(),
                                    residue[atom_chain[-2]].get_vector(),
                                    residue[atom_chain[-1]].get_vector()
                                    )
                bond_len = residue[atom_chain[-1]] - residue[atom_chain[-2]]
                self._atoms[atom_chain[-1]].set_geometry(geom_triple=(dih,
                                                                      planar,
                                                                      bond_len)
                                                         )
            else:
                # if target atom is missing: derive its geometry from siblings
                siblings = self._hierarchy[atom_chain[-2]]
                tgt_sib_pos = [i for i in range(len(siblings)) if siblings[i] == atom_chain[-1]][0]
                sib_dihs = np.full(len(siblings), None)
                sib_planars = np.full(len(siblings), None)
                sib_bond_lens = np.full(len(siblings), None)
                for k in range(len(siblings)):
                    if not residue.has_id(siblings[k]):
                        continue
                    sib_dihs[k] = calc_dihedral(residue[atom_chain[-4]].get_vector(),
                                                residue[atom_chain[-3]].get_vector(),
                                                residue[atom_chain[-2]].get_vector(),
                                                residue[siblings[k]].get_vector()
                                                )
                    sib_planars[k] = calc_angle(residue[atom_chain[-3]].get_vector(),
                                                residue[atom_chain[-2]].get_vector(),
                                                residue[siblings[k]].get_vector()
                                                )
                    sib_bond_lens[k] = residue[siblings[k]] - residue[atom_chain[-2]]
                if np.where(sib_bond_lens == None)[0].shape[0] == len(siblings):
                    raise Exception(f"Atom {atom_chain[-1]} missing is residue. Failed to reproduce geometry!")

                # TODO: fix bond length for different atom pairs?
                bond_len = sib_bond_lens[sib_bond_lens != None].mean()
                planar = sib_planars[sib_planars != None].mean()
                ref_dih_pos = np.where(sib_dihs != None)[0][0]
                dih = sib_dihs[ref_dih_pos] + 2*np.pi/len(siblings) * ((tgt_sib_pos - ref_dih_pos + 3) % 3)
                self._atoms[atom_chain[-1]].set_geometry(geom_triple=(dih,
                                                                      planar,
                                                                      bond_len)
                                                         )
        # setting children
        if atom_chain[-1] not in self._hierarchy.keys():
            return
        children = self._hierarchy[atom_chain[-1]]
        for at_name in children:
            self._init_branch_w_pdb_residue(atom_chain + [at_name], residue)

    def init_dihedrals(self, chi1_value):
        # raise Exception("Not implemented!")
        self.set_chi(1, chi1_value)
        prev_ang_val = chi1_value
        for k in range(1, get_amino_acid_chi_num(self.resname)):
            ang_tgt_atname = get_chi_atoms(self.resname, k + 1)[-1]
            k_ang_pot = get_normalized_rot_potential(self.resname, k + 1, prev_ang_val)
            k_ang_val = np.argmin(k_ang_pot) / k_ang_pot.shape[0] * (2*np.pi)
            prev_ang_val = k_ang_val
            self.set_chi(k + 1, k_ang_val)

    def get_dihedrals(self):
        ans = {}
        # atom_coords = self.get_atom_coords()
        for k in range(get_amino_acid_chi_num(self.resname)):
            # ang = calc_dihedral(*[Vector(*atom_coords[at_name]) for at_name in get_chi_atoms(self.resname, k + 1)])
            ang = self._atoms[get_chi_atoms(self.resname, k + 1)[-1]].get_dihedral()
            ans["chi" + str(k + 1)] = ang
        return ans

    def get_atom_coords(self):
        return self._get_branch_coords("N", np.empty((0,3)))

    def get_volume(self):
        return VolumeProvider().compute_total_volume(self.get_atom_coords())

    def get_hydrogen_status(self):
        return list(self._hydrogen_bonding.keys())

    def _get_branch_coords(self, root_atom_name, prev_atoms):
        res = {}
        root_coord = self._atoms[root_atom_name].get_coord(prev_atoms)
        res[root_atom_name] = root_coord
        if root_atom_name not in self._hierarchy.keys():
            return res
        prev_atoms = np.vstack((prev_atoms, [root_coord]))
        # if len(prev_atoms) > 3:
        prev_atoms = prev_atoms[-3:]
        children = self._hierarchy[root_atom_name]
        for at_name in children:
            res = {**res, **self._get_branch_coords(at_name, prev_atoms)}
        return res

    def set_atom_geometry(self, atom_name, geom_triple=None, fix_coord=None):
        self._atoms[atom_name].set_geometry(geom_triple, fix_coord)

    def _set_atom_dihedral(self, atom_name, dihedral):
        dih_diff = dihedral - self._atoms[atom_name].get_dihedral()
        self._atoms[atom_name].set_dihedral(dihedral)

        # rotating neighbor atom
        parent = [k for k, v in self._hierarchy.items() if atom_name in v][0]
        children = self._hierarchy[parent]
        if len(children) > 1:
            neigh_atom_name = [at for at in children if at != atom_name][0]
            self._atoms[neigh_atom_name].set_dihedral(self._atoms[neigh_atom_name].get_dihedral() + dih_diff)

    def set_chi(self, chi_num, dihedral):
        self._set_atom_dihedral(get_chi_atoms(self.resname, chi_num)[-1], dihedral)
