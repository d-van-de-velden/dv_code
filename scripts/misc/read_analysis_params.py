import json 

def read_analysis_params(fname):

    with open(fname) as f:
        params = json.load(f)

    return params