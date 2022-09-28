from .hydrogen_bonding import HydrogenBonding
import Bio.PDB as pdb
import numpy as np
from .constants import *
from .volume import VolumeProvider
from .hierarchies import hierarchy_dict
import itertools


class IntersectionManager:
    def __init__(self):
        self._hydr_mgr = HydrogenBonding()
        # self._intersect_mgr = IntersectionManager()

    def calc_residue_intersection(self, chain, i,
                                  tgt_res_at_coords,
                                  omit_backbone=True,
                                  omit_sidechain=False,
                                  omit_hydrogens=True,
                                  return_detailed=False,
                                  return_hydrogen_likelihood=False
                                  ):
        substract = 0
        if return_detailed:
            inters_volumes = {}  # key: (f"{j} ({resname})", i_at, j_at), value: intersection volume
        if return_hydrogen_likelihood:
            hydr_bond_pot = 0
        for j in range(len(chain)):
            if i == j:
                continue
            if chain[i][RESIDUE_DIST_REF_ATOM] - chain[j][RESIDUE_DIST_REF_ATOM] > RESIDUE_DIST_THRESHOLD:
                # avoiding too far residues
                continue

            j_hierarchy = hierarchy_dict[chain[j].resname]
            j_sch_atoms = (set(j_hierarchy.keys()) | set(list(itertools.chain(*j_hierarchy.values()))))
            #         j_sch_atoms.remove("CA")
            for i_at in tgt_res_at_coords.keys():
                if i_at in ["N", "CA", "C", "O", "CB"]:
                    continue
                if omit_hydrogens and i_at[0] == "H":
                    continue
                for j_at in chain[j].child_list:
                    if return_hydrogen_likelihood:
                        hydr_bond_pot += self._hydr_mgr.calc_hydrogen_bond_likelihood(i_at,
                                                                       tgt_res_at_coords,
                                                                       chain[i].resname,
                                                                       j_at.name,
                                                                       chain[j])[0]

                    if omit_backbone and (j_at.name in ["N", "CA", "C", "O"]):
                        continue
                    if omit_sidechain and (j_at.name in j_sch_atoms):
                        continue
                    if omit_hydrogens and j_at.name[0] == "H":
                        continue
                    v = VolumeProvider.compute_intersection_volume(tgt_res_at_coords[i_at],
                                                                   j_at.coord,
                                                                   VolumeProvider.vdw_rads[i_at[0]],
                                                                   VolumeProvider.vdw_rads[j_at.name[0]])
                    if v > 0 and return_detailed:
                        inters_volumes[(f"{j} ({chain[j].resname})", i_at, j_at.name)] = v
                    substract += v
        res = [substract]
        if return_detailed:
            #         return substract, inters_volumes
            res += [inters_volumes]
        if return_hydrogen_likelihood:
            res += [hydr_bond_pot]
        if len(res) == 1:
            return res[0]
        else:
            return tuple(res)

