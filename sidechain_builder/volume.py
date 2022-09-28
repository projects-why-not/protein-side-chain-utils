
import numpy as np


class VolumeProvider:
    def __init__(self):
        pass

    vdw_rads = {"H": 1.2,
                "C": 1.7,
                "N": 1.55,
                "O": 1.52,
                "F": 1.47,
                "P": 1.8,
                "S": 1.8
                }

    @classmethod
    def sphere_volume(self, r):
        return 4/3 * np.pi * r**3

    @classmethod
    def compute_intersection_volume(self, o1, o2, r1, r2):
        o1o2 = np.linalg.norm(o2 - o1)
        if o1o2 >= r1 + r2:
            return 0
        if o1o2 <= max(r1, r2) - min(r1, r2):
            # if one sphere is inside another
            return VolumeProvider.sphere_volume(min(r1, r2))

        p = (r1 + r2 + o1o2) / 2
        s_o1o2p = np.sqrt(p * (p - r1) * (p - r2) * (p - o1o2))
        pk = (s_o1o2p * 2 / o1o2)
        h1 = r1 - np.sqrt(r1 ** 2 - pk ** 2)
        h2 = r2 - np.sqrt(r2 ** 2 - pk ** 2)
        v = 2 / 3 * np.pi * r1 ** 2 * h1 + 2 / 3 * np.pi * r2 ** 2 * h2
        v -= 1 / 3 * np.pi * pk ** 2 * o1o2
        return v

    @classmethod
    def compute_total_volume(self, atom_coords):
        v = 0
        for at_name in atom_coords:
            v += VolumeProvider.sphere_volume(self.vdw_rads[at_name[0]])

        at_names = list(atom_coords.keys())
        for i in range(len(at_names) - 1):
            for j in range(i+1, len(at_names)):
                v -= VolumeProvider.compute_intersection_volume(atom_coords[at_names[i]],
                                                      atom_coords[at_names[j]],
                                                      VolumeProvider.vdw_rads[at_names[i][0]],
                                                      VolumeProvider.vdw_rads[at_names[j][0]]
                                                      )
        return v
