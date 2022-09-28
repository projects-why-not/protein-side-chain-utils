import numpy as np
import os
import json

ROTAMER_POTENTIAL_GRID_DEGREE_STEP = 20


def get_normalized_rot_potential(resname, chi_num, prev_ang_value):
    """
    :param resname: name of amino acid.
    :param chi_num: number of chi angle (starting from 1).
    :param prev_ang_value: value for previous angle along the chain ([0, 2*pi)). If chi_num is 1, this
                            should be a pair (phi, psi), phi, psi in [-pi, pi].
    :return: normalized potential for chi_{chi_num}.
    """
    if chi_num > 1:
        prev_ang_value = (prev_ang_value + 2*np.pi) % (2*np.pi)
    else:
        prev_ang_value = np.array(prev_ang_value)

    data = np.load("/".join(os.path.realpath(__file__).split("/")[:-2]) +
                   "/resources/rotameric_potentials/{}.npy".format(resname),
                   allow_pickle=True)
    ang_data = data.item()["chi_" + str(chi_num)]
    # print(ang_data.shape)
    # input()
    step = ROTAMER_POTENTIAL_GRID_DEGREE_STEP / 180 * np.pi
    if chi_num == 1:
        if prev_ang_value[0] is None:
            pot = ang_data[:, int(prev_ang_value[1] // step)].mean(axis=0)
            return (pot - pot.min()) / (pot.max() - pot.min())
        elif prev_ang_value[1] is None:
            pot = ang_data[int(prev_ang_value[0] // step)].mean(axis=0)
            return (pot - pot.min()) / (pot.max() - pot.min())

        coord = (prev_ang_value // step).astype(int)
        # coord = (prev_ang_value // step).astype(int)
        # coord += ang_data.shape[0] // 2
        # print("\n", prev_ang_value)
        # print(coord)
        # coord = np.maximum(np.minimum(coord, ang_data.shape[0] - 1), 0)
        # print(coord)
        # input()
        return ang_data[coord[0], coord[1]]
    else:
        coord = int(prev_ang_value // step)
        return ang_data[coord]

    # res = None
    # with open("/".join(os.path.realpath(__file__).split("/")[:-2]) + "/resources/normed_rot_potentials.json") as f:
    #     res = json.loads(f.read())[resname]
    # return np.array(res)


if __name__ == "__main__":
    print(get_normalized_rot_potential("VAL", 1, (-np.pi, np.pi)))
    # print(get_normalized_rot_potential("PRO", 2, 20/180*np.pi))
