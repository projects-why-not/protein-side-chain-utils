import json
import os


def get_agg_pot_params(resname, kind="w_all"):
    res = None
    with open("/".join(os.path.realpath(__file__).split("/")[:-2]) + "/resources/agg_pot_params.json") as f:
        res = json.loads(f.read())[resname][kind]
    return res
