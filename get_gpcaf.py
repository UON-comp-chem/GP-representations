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
from gplearn import coordination_features
from gplearn.defaults import elements_list, elements_number
from copy import deepcopy
import random

db_path = './data/scidata-original-oxygen0328.json'
with open(db_path, 'r') as file_handle:
     db = json.load(file_handle)
db = db[:10]
y = np.array([n['energy'] for n in db])
a = np.arange(len(db))
np.random.shuffle(a)

inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate='O')
outer_fingerprinter = coordination_features.OuterShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate='O')
caf = coordination_features.StackedFingerprinter(inner_fingerprinter, outer_fingerprinter)

X = caf.fit_transform(db)
y = np.array([n['energy'] for n in db])
X = X.astype(np.float32)

print("X:\n", X)
print("y:\n", y)