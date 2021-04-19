#!/usr/bin/env python
# coding: utf-8
# Author Xinyu Li

import sys
import os
# Insert system path to use scripts in gplearn
sys.path.append(os.path.dirname(os.getcwd()))

import json
import numpy as np

from gplearn.read_datasets import make_initial_atoms_from_doc, make_atoms_from_doc
from gplearn.gpfchl import GPFCHLFingerprinter
from gplearn.defaults import elements_list, elements_number

from copy import deepcopy
import random

rdb_path = './data/scidata-rebuild-oxygen0328.json' # regression data set
params = {"nRs2": 21, 
            "nRs3": 26, 
            "eta2": 0.21802291630597279, 
            "eta3": 1.5053496197211147, 
            "two_body_decay": 2.078515163355548, 
            "three_body_decay": 1.9426424314522792, 
            "three_body_weight": 36.0107794486399,
            "ads_idx":12}

with open(rdb_path, 'r') as file_handle:
    rdb = json.load(file_handle)
rdb = rdb[:10]
fchl = GPFCHLFingerprinter(elements= elements_number, relaxed = True, rcut = 6.0, acut = 6.0, ads_idx = 12, **params)

y = np.array([n['energy'] for n in rdb])
X = fchl.fit_transform(rdb, y)
X = X.astype(np.float32)

print("X:\n", X)
print("y:\n", y)