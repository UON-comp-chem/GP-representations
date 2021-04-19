import os
import sys
sys.path.append(os.path.dirname(os.getcwd()))


import json
import numpy as np

from gplearn.read_datasets import make_initial_atoms_from_doc, make_atoms_from_doc
from gplearn import coordination_features
from gplearn.defaults import elements_list, elements_number
from copy import deepcopy
import random


    
def test_gp():
    docs_gen = []
    for ele in elements_list:
        _d = {"slab_name":ele,
               "metalA":ele,
               "metalB":"nan",
               "coordination":ele,
               "nextnearest":ele}
        docs_gen.append(_d)
    inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["group_id","period"],elements = elements_list, adsorbate='O')
    outer_fingerprinter = coordination_features.NextNearestFingerprinter(features = ["group_id","period"],elements = elements_list, adsorbate='O')
    caf = coordination_features.StackedFingerprinter(inner_fingerprinter, outer_fingerprinter)

    X = caf.fit_transform(docs_gen)
    from ase.data import atomic_numbers
    from gplearn.defaults import PTP
    import mendeleev
    
    for idx, ele in enumerate(elements_list):
        gp = X[idx][:2]
        assert gp[0] == getattr(mendeleev, ele).group_id
        assert gp[1] == getattr(mendeleev, ele).period
        

def test_an():
    docs_gen = []
    for ele in elements_list:
        _d = {"slab_name":ele,
               "metalA":ele,
               "metalB":"nan",
               "coordination":ele,
               "nextnearest":ele}
        docs_gen.append(_d)
    inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["atomic_number"],elements = elements_list, adsorbate='O')
    outer_fingerprinter = coordination_features.NextNearestFingerprinter(features = ["atomic_number"],elements = elements_list, adsorbate='O')
    caf = coordination_features.StackedFingerprinter(inner_fingerprinter, outer_fingerprinter)

    X = caf.fit_transform(docs_gen)
    from ase.data import atomic_numbers
    from gplearn.defaults import PTP
    import mendeleev
    
    for idx, ele in enumerate(elements_list):
        print(ele)
        assert X[idx][0] == atomic_numbers[ele]
        

    
def test_count():
    ele = "Pt"
    docs_gen = []
    for c in range(1, 20):
        _d = {"slab_name":ele,
               "metalA":ele,
               "metalB":"nan",
               "coordination":((ele+"-") * c)[:-1],
               "nextnearest":((ele+"-") * c)[:-1]}
        
        docs_gen.append(_d)
    inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["count"],elements = elements_list, adsorbate='O')
    outer_fingerprinter = coordination_features.NextNearestFingerprinter(features = ["count"],elements = elements_list, adsorbate='O')
    caf = coordination_features.StackedFingerprinter(inner_fingerprinter, outer_fingerprinter)

    X = caf.fit_transform(docs_gen)
    from ase.data import atomic_numbers
    from gplearn.defaults import PTP
    import mendeleev
    
    for idx, c in enumerate(range(1, 20)):
        assert X[idx][0] == c
                                                                                     
                                                                                     