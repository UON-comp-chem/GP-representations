#!/usr/bin/env python
# coding: utf-8
# Author Xinyu Li

import sys
import os
# Insert system path to use scripts in gplearn
sys.path.append(os.path.dirname(os.getcwd()))

import json
import pickle
import numpy as np

from gplearn.read_datasets import make_initial_atoms_from_doc, make_atoms_from_doc
from gplearn import coordination_features
from gplearn.defaults import elements_list, elements_number
from gplearn.atoms_operators import fingerprint_adslab
from copy import deepcopy
import random
from tqdm import tqdm

db_path = './data/scidata-hydrogen-rebuilt.pkl'
ADSORBATE = "H"

with open(db_path, 'rb') as file_handle:
     db = pickle.load(file_handle)
for doc in tqdm(db):
    atoms = make_atoms_from_doc(doc)
    symbols = [atom.symbol for atom in atoms[atoms.get_tags() == 0]]
    unique = np.unique(symbols)
    slab_name  = "".join(unique)
    if len(unique) == 2:
        metalA, metalB= unique[0], unique[1]
    elif len(unique) == 1:
        metalA, metalB = unique[0], 'nan'
    else:
        raise "only bimetallic alloys supported, please check"
    doc['slab_name'] = slab_name
    doc['metalA'] = metalA
    doc['metalB'] = metalB
    ffp = fingerprint_adslab(atoms, adsorbate = ADSORBATE[0])
    for name, value in ffp.items():
        doc[name] = value

db = db[:10]
y = np.array([n['energy'] for n in db])
a = np.arange(len(db))
np.random.shuffle(a)

# inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate='O')
# outer_fingerprinter = coordination_features.OuterShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate='O')

inner_fingerprinter = coordination_features.InnerShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate=ADSORBATE[0])
outer_fingerprinter = coordination_features.OuterShellFingerprinter(features = ["group_id","period", "median_energy", "electronegativity", "count"],elements = elements_list, adsorbate=ADSORBATE[0])

caf = coordination_features.StackedFingerprinter(inner_fingerprinter, outer_fingerprinter)

X = caf.fit_transform(db)
y = np.array([n['energy'] for n in db])
X = X.astype(np.float32)

print("X:\n", X)
print("y:\n", y)