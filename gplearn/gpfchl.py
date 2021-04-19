#!/usr/bin/env python
# coding: utf-8
# MIT License

"""
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.

---
This code borrows heavily from 
https://github.com/qmlcode/qml/blob/develop/qml/qmlearn/representations.py

License: 
---

#
# Copyright (c) 2018 Lars Andersen Bratholm
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""


from .mongo import make_initial_atoms_from_doc, make_atoms_from_doc
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from qml.utils.alchemy import PTP
from .fgpfchl.representations import (generate_local_two_body_fchl_acsf, 
                                  generate_local_gp_fchl, 
                                  generate_local_fchl_acsf,
                                  generate_local_two_body_gp_fchl)

def _get_supercell(obj, rCut=5.0):
    
    rCutHard = rCut #+ 5  # Giving extra space for hard cutOff
    """ 
    Shameless steal from: https://github.com/SINGROUP/SOAPLite/blob/master/soaplite/core.py
    Credit goes to SOAPlite developers
    Takes atoms object (with a defined cell) and a radial cutoff.
    Returns a supercell centered around the original cell
    generously extended to contain all spheres with the given radial
    cutoff centered around the original atoms.
    Args:
        obj                   The `ase.Atoms` object with PBC that you want
                              to get its supercell
    Returns:
        shifted_suce          A `ase.Atoms` object which is obj's supercell, 
                              it preserves the chemical environment of all atoms
                              in the obj
    """
    cell_vectors = obj.get_cell()
    a1, a2, a3 = cell_vectors[0], cell_vectors[1], cell_vectors[2]

    # vectors perpendicular to two cell vectors
    b1 = np.cross(a2, a3, axis=0)
    b2 = np.cross(a3, a1, axis=0)
    b3 = np.cross(a1, a2, axis=0)
    # projections onto perpendicular vectors
    p1 = np.dot(a1, b1) / np.dot(b1, b1) * b1
    p2 = np.dot(a2, b2) / np.dot(b2, b2) * b2
    p3 = np.dot(a3, b3) / np.dot(b3, b3) * b3
    xyz_arr = np.linalg.norm(np.array([p1, p2, p3]), axis=1)
    cell_images = np.ceil(rCutHard/xyz_arr)
    nx = int(cell_images[0])
    ny = int(cell_images[1])
    #nz = int(cell_images[2])
    suce = obj * (1+2*nx, 1+2*ny, 1)
    shift = obj.get_cell()

    shifted_suce = suce.copy()
    shifted_suce.translate(-shift[0]*nx - shift[1]*ny)
    start = int(len(obj) *((1+2*ny)*nx + ny))
    end = int(start + len(obj))
    indices = np.concatenate((np.arange(start, end), np.arange(0, start), np.arange(end, len(obj) *(1+2*ny)*(1 + 2*nx))))
    shifted_suce = shifted_suce[indices]
    return shifted_suce


class FCHLFingerprinter(BaseEstimator, TransformerMixin):
    """
    The variant of Atom-Centered Symmetry Functions 
    :param data: Optional Data object containing all molecules used in training \
            and/or prediction
    :type data: Data object
    :param elements: Atomnumber of elements that the representation should support.
                     `elements='auto'` will try to determine this automatically.
    :type elements: list
    """

    def __init__(self, elements='auto',relaxed = False,
            nRs2=24, nRs3=20, nFourier=1, eta2=0.32, eta3=2.7, zeta=np.pi, rcut=6.0, acut=6.0,
            two_body_decay=1.8, three_body_decay=0.57, three_body_weight=13.4, ads_idx = -1):
        super(FCHLFingerprinter,self).__init__()
        self.nRs2 = nRs2
        self.nRs3 = nRs3
        self.nFourier = nFourier
        self.eta2 = eta2 
        self.eta3 = eta3
        self.zeta = zeta
        self.rcut = rcut
        self.acut = acut
        self.two_body_decay = two_body_decay
        self.three_body_decay = three_body_decay
        self.three_body_weight = three_body_weight
        self.elements = elements
        self.relaxed = relaxed
        self.ads_idx = ads_idx

    def fit(self, docs, y):
        if self.relaxed == False:
            self.structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            self.structures = [make_atoms_from_doc(n) for n in docs]
        fitted_ele = np.unique([m for n in self.structures for m in n.get_atomic_numbers()])
        if isinstance(self.elements, str):
            if self.elements == 'auto':
                self.elements = fitted_ele
        else:
            diff = np.setdiff1d(fitted_ele, self.elements) 
            if len(diff) != 0: # If some elements not included, then add it 
                self.elements = np.concatenate((self.elements, diff))
        return self

    def transform(self, docs):
        if self.relaxed == False:
            structures = [make_initial_atoms_from_doc(n) for n in docs]
        elif self.relaxed == True:
            structures = [make_atoms_from_doc(n) for n in docs]
        natoms = [len(atoms) for atoms in structures]
        
        structures = [_get_supercell(obj, rCut=self.rcut) for obj in structures]
        if self.ads_idx >= 0:
            ads_idx = [self.ads_idx for n in natoms]
        else:
            ads_idx = [n + self.ads_idx for n in natoms]
        
        # quike way on scientific data
        #structures = [atoms *(5, 5, 1) for atoms in structures]
        #if self.ads_idx >= 0:
        #    ads_idx = [int(12 * n + self.ads_idx) for n in natoms]
        #else:
        #    ads_idx = [int(13 * n + self.ads_idx) for n in natoms]

        assert len(np.unique([atoms[idx].symbol for atoms, idx in zip(structures, ads_idx)])) == 1
        nuclear_charges = np.array([m.get_atomic_numbers() for m in structures])
        coordinates = np.array([m.positions for m in structures])
        natoms = np.array([len(m) for m in nuclear_charges])
        max_atoms = np.amax(natoms)

        self._check_elements(nuclear_charges)

        repo = generate_local_fchl_acsf(ads_idx[0], nuclear_charges[0], coordinates[0] , elements=self.elements,
                           nRs2=self.nRs2, nRs3=self.nRs3, nFourier=self.nFourier, eta2=self.eta2, 
                           eta3=self.eta3, zeta=self.zeta, rcut=self.rcut, acut=self.acut,
                           two_body_decay=self.two_body_decay, three_body_decay=self.three_body_decay, 
                           three_body_weight=self.three_body_weight)
        representations = np.zeros((len(natoms), len(repo)))
        i = 0
        for idx, charge, xyz in zip(ads_idx, nuclear_charges, coordinates):
            rep = generate_local_fchl_acsf(idx, charge, xyz, elements=self.elements,
                           nRs2=self.nRs2, nRs3=self.nRs3, nFourier=self.nFourier, eta2=self.eta2, 
                           eta3=self.eta3, zeta=self.zeta, rcut=self.rcut, acut=self.acut,
                           two_body_decay=self.two_body_decay, three_body_decay=self.three_body_decay, 
                           three_body_weight=self.three_body_weight)
            representations[i] = rep
            i += 1
        return representations
    
    def _check_elements(self, nuclear_charges):
        elements_transform = np.unique([m for n in nuclear_charges for m in n])
        if not np.isin(elements_transform, self.elements).all():
            print("Warning: Trying to transform molecules with elements",
                  "not included during fit in the %s method." % self.__class__.__name__,
                  "%s used in training but trying to transform %s" % (str(self.elements), str(elements_transform)))

class TwoBodyFCHLFingerprinter(BaseEstimator, TransformerMixin):
    """
    The variant of Atom-Centered Symmetry Functions 
    :param data: Optional Data object containing all molecules used in training \
            and/or prediction
    :type data: Data object
    :param elements: Atomnumber of elements that the representation should support.
                     `elements='auto'` will try to determine this automatically.
    :type elements: list
    """

    def __init__(self, elements='auto',relaxed = False,
            nRs2=24, eta2=0.32,  rcut=6.0, two_body_decay=1.8, ads_idx = -1):
        super(TwoBodyFCHLFingerprinter,self).__init__()
        self.nRs2 = nRs2
        self.eta2 = eta2 
        self.rcut = rcut
        self.two_body_decay = two_body_decay
        self.elements = elements
        self.relaxed = relaxed
        self.ads_idx = ads_idx

    def fit(self, docs, y):
        if self.relaxed == False:
            self.structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            self.structures = [make_atoms_from_doc(n) for n in docs]
        fitted_ele = np.unique([m for n in self.structures for m in n.get_atomic_numbers()])
        if isinstance(self.elements, str):
            if self.elements == 'auto':
                self.elements = fitted_ele
        else:
            diff = np.setdiff1d(fitted_ele, self.elements)
            if len(diff) != 0:
                self.elements = np.concatenate((self.elements, diff))
        return self

    def transform(self, docs):
        if self.relaxed == False:
            structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            structures = [make_atoms_from_doc(n) for n in docs]
        natoms = [len(atoms) for atoms in structures]
        
        structures = [_get_supercell(obj, rCut=self.rcut) for obj in structures]
        if self.ads_idx >= 0:
            ads_idx = [self.ads_idx for n in natoms]
        else:
            ads_idx = [n + self.ads_idx for n in natoms]
        
        # structures = [atoms *(5, 5, 1) for atoms in structures]
        # if self.ads_idx >= 0:
            # ads_idx = [int(12 * n + self.ads_idx) for n in natoms]
        # else:
            # ads_idx = [int(13 * n + self.ads_idx) for n in natoms]
        
        assert len(np.unique([atoms[idx].symbol for atoms, idx in zip(structures, ads_idx)])) == 1
        nuclear_charges = np.array([m.get_atomic_numbers() for m in structures])
        coordinates = np.array([m.positions for m in structures])
        natoms = np.array([len(m) for m in nuclear_charges])
        max_atoms = np.amax(natoms)

        self._check_elements(nuclear_charges)

        repo = generate_local_two_body_fchl_acsf(ads_idx[0], nuclear_charges[0], coordinates[0], elements = self.elements,
                        nRs2=self.nRs2,  eta2=self.eta2, rcut=self.rcut, two_body_decay=self.two_body_decay, pad=False)
        representations = np.zeros((len(natoms), len(repo)))
        i = 0
        for idx, charge, xyz in zip(ads_idx, nuclear_charges, coordinates):
            rep = generate_local_two_body_fchl_acsf(idx, charge, xyz, elements = self.elements,
                        nRs2=self.nRs2,  eta2=self.eta2, rcut=self.rcut, two_body_decay=self.two_body_decay, pad=False)
            representations[i] = rep
            i += 1
        return representations
    
    def _check_elements(self, nuclear_charges):
        elements_transform = np.unique([m for n in nuclear_charges for m in n])
        if not np.isin(elements_transform, self.elements).all():
            print("Warning: Trying to transform molecules with elements",
                  "not included during fit in the %s method." % self.__class__.__name__,
                  "%s used in training but trying to transform %s" % (str(self.elements), str(elements_transform)))

class GPFCHLFingerprinter(BaseEstimator, TransformerMixin):
    """
    The variant of Atom-Centered Symmetry Functions 
    :param data: Optional Data object containing all molecules used in training \
            and/or prediction
    :type data: Data object
    :param elements: Atomnumber of elements that the representation should support.
                     `elements='auto'` will try to determine this automatically.
    :type elements: list
    """

    def __init__(self, elements='auto',relaxed = False,
            nRs2=24, nRs3=20, nFourier=1, eta2=0.32, eta3=2.7, zeta=np.pi, rcut=6.0, acut=6.0,
            two_body_decay=1.8, three_body_decay=0.57, three_body_weight=13.4, ads_idx = -1):
        super(GPFCHLFingerprinter,self).__init__()
        self.nRs2 = nRs2
        self.nRs3 = nRs3
        self.nFourier = nFourier
        self.eta2 = eta2 
        self.eta3 = eta3
        self.zeta = zeta
        self.rcut = rcut
        self.acut = acut
        self.two_body_decay = two_body_decay
        self.three_body_decay = three_body_decay
        self.three_body_weight = three_body_weight
        self.elements = elements
        self.relaxed = relaxed
        self.ads_idx = ads_idx

    def fit(self, docs, y):
        if self.relaxed == False:
            self.structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            self.structures = [make_atoms_from_doc(n) for n in docs]
        fitted_ele = np.unique([m for n in self.structures for m in n.get_atomic_numbers()])
        if isinstance(self.elements, str):
            if self.elements == 'auto':
                self.elements = fitted_ele
        else:
            diff = np.setdiff1d(fitted_ele, self.elements)
            if len(diff) != 0:
                self.elements = np.concatenate((self.elements, diff))

        return self

    def transform(self, docs):
        if self.relaxed == False:
            structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            structures = [make_atoms_from_doc(n) for n in docs]
        natoms = [len(atoms) for atoms in structures]
        
        structures = [_get_supercell(obj, rCut=self.rcut) for obj in structures]
        if self.ads_idx >= 0:
            ads_idx = [self.ads_idx for n in natoms]
        else:
            ads_idx = [n + self.ads_idx for n in natoms]
        
        # structures = [atoms *(5, 5, 1) for atoms in structures]
        # if self.ads_idx >= 0:
            # ads_idx = [int(12 * n + self.ads_idx) for n in natoms]
        # else:
            # ads_idx = [int(13 * n + self.ads_idx) for n in natoms]
        assert len(np.unique([atoms[idx].symbol for atoms, idx in zip(structures, ads_idx)])) == 1
        nuclear_charges = np.array([m.get_atomic_numbers() for m in structures])
        coordinates = np.array([m.positions for m in structures])
        natoms = np.array([len(m) for m in nuclear_charges])
        max_atoms = np.amax(natoms)

        self._check_elements(nuclear_charges)

        repo = generate_local_gp_fchl(ads_idx[0], nuclear_charges[0], coordinates[0] , elements=self.elements,
                           nRs2=self.nRs2, nRs3=self.nRs3, nFourier=self.nFourier, eta2=self.eta2, 
                           eta3=self.eta3, zeta=self.zeta, rcut=self.rcut, acut=self.acut,
                           two_body_decay=self.two_body_decay, three_body_decay=self.three_body_decay, 
                           three_body_weight=self.three_body_weight,
                           pad=False)
        representations = np.zeros((len(natoms), len(repo)))
        i = 0
        for idx, charge, xyz in zip(ads_idx, nuclear_charges, coordinates):
            rep = generate_local_gp_fchl(idx, charge, xyz , elements=self.elements,
                           nRs2=self.nRs2, nRs3=self.nRs3, nFourier=self.nFourier, eta2=self.eta2, 
                           eta3=self.eta3, zeta=self.zeta, rcut=self.rcut, acut=self.acut,
                           two_body_decay=self.two_body_decay, three_body_decay=self.three_body_decay, 
                           three_body_weight=self.three_body_weight,
                           pad=False)
            representations[i] = rep
            i += 1
        return representations
    
    def _check_elements(self, nuclear_charges):
        elements_transform = np.unique([m for n in nuclear_charges for m in n])
        if not np.isin(elements_transform, self.elements).all():
            print("Warning: Trying to transform molecules with elements",
                  "not included during fit in the %s method." % self.__class__.__name__,
                  "%s used in training but trying to transform %s" % (str(self.elements), str(elements_transform)))
            
class TwoBodyGPFCHLFingerprinter(BaseEstimator, TransformerMixin):
    """
    The variant of Atom-Centered Symmetry Functions 
    :param data: Optional Data object containing all molecules used in training \
            and/or prediction
    :type data: Data object
    :param elements: Atomnumber of elements that the representation should support.
                     `elements='auto'` will try to determine this automatically.
    :type elements: list
    """

    def __init__(self, elements='auto',relaxed = False,
            nRs2=24, eta2=0.32, rcut=6.0, two_body_decay=1.8, ads_idx = -1):
        super(TwoBodyGPFCHLFingerprinter,self).__init__()
        self.nRs2 = nRs2
        self.eta2 = eta2
        self.rcut = rcut
        self.two_body_decay = two_body_decay
        self.elements = elements
        self.relaxed = relaxed
        self.ads_idx = ads_idx

    def fit(self, docs, y):
        if self.relaxed == False:
            self.structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            self.structures = [make_atoms_from_doc(n) for n in docs]
        fitted_ele = np.unique([m for n in self.structures for m in n.get_atomic_numbers()])
        if isinstance(self.elements, str):
            if self.elements == 'auto':
                self.elements = fitted_ele
        else:
            diff = np.setdiff1d(fitted_ele, self.elements)
            if len(diff) != 0:
                self.elements = np.concatenate((self.elements, diff))

        return self

    def transform(self, docs):
        if self.relaxed == False:
            structures = [make_initial_atoms_from_doc(n) for n in docs]
        else:
            structures = [make_atoms_from_doc(n) for n in docs]
        natoms = [len(atoms) for atoms in structures]
        
        structures = [_get_supercell(obj, rCut=self.rcut) for obj in structures]
        if self.ads_idx >= 0:
            ads_idx = [self.ads_idx for n in natoms]
        else:
            ads_idx = [n + self.ads_idx for n in natoms]
        
        # structures = [atoms *(5, 5, 1) for atoms in structures]
        # if self.ads_idx >= 0:
            # ads_idx = [int(12 * n + self.ads_idx) for n in natoms]
        # else:
            # ads_idx = [int(13 * n + self.ads_idx) for n in natoms]
        
        assert len(np.unique([atoms[idx].symbol for atoms, idx in zip(structures, ads_idx)])) == 1
        nuclear_charges = np.array([m.get_atomic_numbers() for m in structures])
        coordinates = np.array([m.positions for m in structures])
        natoms = np.array([len(m) for m in nuclear_charges])
        max_atoms = np.amax(natoms)

        self._check_elements(nuclear_charges)

        repo = generate_local_two_body_gp_fchl(ads_idx[0], nuclear_charges[0], coordinates[0] , elements=self.elements,
                           nRs2=self.nRs2, eta2=self.eta2, rcut=self.rcut, 
                           two_body_decay=self.two_body_decay, pad=False)
        representations = np.zeros((len(natoms), len(repo)))
        i = 0
        for idx, charge, xyz in zip(ads_idx, nuclear_charges, coordinates):
            rep = generate_local_two_body_gp_fchl(idx, charge, xyz, elements=self.elements,
                           nRs2=self.nRs2, eta2=self.eta2, rcut=self.rcut, 
                           two_body_decay=self.two_body_decay, pad=False)
            representations[i] = rep
            i += 1
        return representations
    
    def _check_elements(self, nuclear_charges):
        elements_transform = np.unique([m for n in nuclear_charges for m in n])
        if not np.isin(elements_transform, self.elements).all():
            print("Warning: Trying to transform molecules with elements",
                  "not included during fit in the %s method." % self.__class__.__name__,
                  "%s used in training but trying to transform %s" % (str(self.elements), str(elements_transform)))
            


