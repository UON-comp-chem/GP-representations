#!/usr/bin/env python
# coding: utf-8
# GNU Lesser General Public License v3.0

"""
This source code is licensed under the GNU Lesser General Public License v3.0
found in the LICENSE file in the root directory of this source tree.

---
This code simply copy from 
https://github.com/ulissigroup/GASpy/blob/master/gaspy/mongo.py
We add it here to solve dependency issues.

All credit goes to the original authors: 
'John Kitchin', 'Kevin Tran'
"""

'''
Note that this entire module was simply copied from John Kitchin's
vasp.mongo for use here. We used his code and put it here to solve
dependency issues. All credit goes to him, and we thank him for his help.
This module will be like the ase-db but different in the following ways:
1. Booleans are stored as booleans.
2. There is no numeric id.
3. Tags are stored in an array.
'''

__authors__ = ['John Kitchin', 'Kevin Tran']
__email__ = 'ktran@andrew.cmu.edu'

from ase.data import chemical_symbols
from ase import Atoms, Atom
from ase.visualize import view
from ase.constraints import dict2constraint
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.data import chemical_symbols
import numpy as np
import os
import csv
from collections import OrderedDict, defaultdict
import datetime
import json
import spglib
import pickle
from tqdm import tqdm
from ase.io.jsonio import encode
from ase.constraints import dict2constraint
from .atoms_operators import fingerprint_adslab
from .coordination_features import CACHE_LOCATION
from .defaults import elements_number, elements_list
from .util import CovalentRadius, pbcs
from scipy.spatial.distance import cdist
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.io import write
from .mongo import make_atoms_from_doc

def remove_trimetallic(db, adsorbate = 'H'):
    if adsorbate == 'H':
        adsnc = [1]
    elif adsorbate == 'CO':
        adsnc = [6, 8]
    newdoc = []
    for doc in db:
        atoms = make_atoms_from_doc(doc)
        an = atoms.get_atomic_numbers()
        an = np.unique(an)
        an = np.setdiff1d(an, adsnc)
        if len(an) <= 2:
            newdoc.append(doc)
    return newdoc

def get_metallics(db, adsorbate = 'CO'):
    path = os.path.join(CACHE_LOCATION, "mp_comp_data.pkl")
    with open(path, 'rb') as fp:
        mpids = pickle.load(fp)
    keepmpids = []
    for key, value in mpids.items():
        diff =  np.setdiff1d(value, elements_list)
        if len(diff) == 0:
            keepmpids.append(key)
    newdoc = []
    if adsorbate == 'CO':
        elements = np.concatenate([elements_number, [6, 8]])
    elif adsorbate == 'H':
        elements = np.concatenate([elements_number, [1]])
    for doc in db:
        if doc['mpid'] in keepmpids:
            atoms = make_atoms_from_doc(doc)
            an = np.unique(atoms.get_atomic_numbers())
            diff = np.setdiff1d(an, elements)
            if len(diff) == 0:
                newdoc.append(doc)
    return newdoc 

def remove_too_much_samples(db, line = 5):
    with open("../data/mp_comp_data.pkl", 'rb') as fp:
        mpid2ele = pickle.load(fp)
    atoms  = [make_atoms_from_doc(doc) for doc in db]
    elements = np.array([chemical_symbols[n] for n in np.unique([n for atom in atoms for n in atom.get_atomic_numbers()])])
    no = np.zeros(len(elements))
    for doc in db:
        _ele = mpid2ele[doc['mpid']]
        for __e in _ele:
            no[elements == __e] += 1
    print(elements, no)
    elements_2b_reduced = [n for n in elements[np.argsort(no)[::-1]][:line]]
    keep_number = no[np.argsort(no)[::-1]].astype(int)
    keep_number = keep_number[line]
    print(elements_2b_reduced, keep_number)
    for ele in elements_2b_reduced:
        ndoc = []
        edoc = []
        for doc in db:
            _ele = mpid2ele[doc['mpid']]
            if ele not in _ele:
                ndoc.append(doc)
            else:
                edoc.append(doc)
        a = np.arange(0, len(edoc))
        np.random.shuffle(a)
        for i in a[:keep_number]:
            ndoc.append(edoc[i])
        db = ndoc
    return db

mat = StructureMatcher(primitive_cell=False)

def remove_non_interaction(data,adsorbate = 'H',upperbound = 1.15,radius=12.0,nnbr=9):
    if adsorbate == 'H':
        adsan = 1
    elif adsorbate == 'CO':
        adsan = 6
    new_data = []
    for datum in tqdm(data):
        if datum['adsorption_energy'] < -3.0 or datum['adsorption_energy'] > 1.0:
            continue
        ## Find binding site atoms
        ### Get atoms
        ai = make_atoms_from_doc(datum['initial_configuration'])
        a  = make_atoms_from_doc(datum)
        ### Adsorbate index
        adsidx = np.where(a.get_tags()==1)[0]
        ### surface index
        surfidx = np.where(a.get_tags()==0)[0]
        surfmap = {i:surfidx[i] for i in range(len(surfidx))}
        an = a.get_atomic_numbers()
        sym = a.get_chemical_symbols()
        ### carbon
        Cidx = adsidx[an[adsidx]==adsan][0]
        ### positions
        CPos = a.get_positions()[Cidx,:]
        SurfPos = a.get_scaled_positions(wrap=False)[surfidx]
        SurfPosPbcs = np.concatenate(np.dot(SurfPos[None,:,:]+np.array(pbcs)[:,None,:],a.cell))
        ### find dist
        dists = cdist([CPos],SurfPosPbcs)[0]
        ### get criteria
        criteria =[]
        for _ in range(len(pbcs)):
            for i in surfmap:
                criteria.append(CovalentRadius[sym[surfmap[i]]]+CovalentRadius['C'])
        criteria = np.array(criteria)*upperbound
        ### apply criteria
        BindingSiteAtomIdx = [surfmap[i] for i in np.remainder(np.where(dists<criteria)[0],len(surfmap))] # map back to original
        if len(BindingSiteAtomIdx) == 0:
            continue # no binding site.
        new_data.append(datum)
    return new_data