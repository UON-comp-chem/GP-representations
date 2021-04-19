#!/usr/bin/env python
# coding: utf-8
# MIT License

"""
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.

---
Shameless steal from: 
https://github.com/qmlcode/qml/blob/develop/qml/representations/representations.py
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


import numpy as np
import itertools as itl
from ..defaults import PTP, group
from .frepresentations import (fgenerate_local_fchl_acsf, 
                              fgenerate_local_two_body_fchl_acsf, 
                              fgenerate_local_gp_fchl, 
                              fgenerate_local_two_body_gp_fchl)

def generate_local_gp_fchl(idx, nuclear_charges, coordinates, elements = [1,6,7,8,16],
                             nRs2=24, nRs3=20, nFourier=1, eta2=0.32, eta3=2.7, zeta=np.pi, rcut=6.0, acut=6.0,
                             two_body_decay=1.8, three_body_decay=0.57, three_body_weight=13.4,
                             pad=False):

    Rs2 = np.linspace(0, rcut, 1+nRs2)[1:]
    Rs3 = np.linspace(0, acut, 1+nRs3)[1:]
    unique_groups = np.sort(np.unique([group(n) for n in elements]))
    if idx < 0:
        idx += len(nuclear_charges)
    Ts = np.linspace(0, np.pi, 2*nFourier)
    n_periods = len(unique_groups)
    natoms = len(coordinates)
    groups = [PTP[n][1] for n in nuclear_charges]
    periods = [PTP[n][0] for n in nuclear_charges]
    

    descr_size = n_periods * nRs2 + (n_periods * (n_periods + 1)) * nRs3* nFourier

    # Normalization constant for three-body 
    three_body_weight = np.sqrt(eta3/np.pi) * three_body_weight

    rep = fgenerate_local_gp_fchl(idx, coordinates, nuclear_charges, periods, groups, unique_groups, Rs2, Rs3,
            Ts, eta2, eta3, zeta, rcut, acut, natoms, descr_size,
            two_body_decay, three_body_decay, three_body_weight)

    return rep

def generate_local_two_body_gp_fchl(idx, nuclear_charges, coordinates, elements = [1,6,7,8,16],
                             nRs2=24, eta2=0.32, rcut=6.0, two_body_decay=1.8, pad=False):

    Rs2 = np.linspace(0, rcut, 1+nRs2)[1:]
    unique_groups = np.sort(np.unique([group(n) for n in elements]))
    if idx < 0:
        idx += len(nuclear_charges)
    n_groups = len(unique_groups)
    natoms = len(coordinates)
    groups = [PTP[n][1] for n in nuclear_charges]
    periods = [PTP[n][0] for n in nuclear_charges]

    descr_size = n_groups * nRs2 

    rep = fgenerate_local_two_body_gp_fchl(idx, coordinates, nuclear_charges, periods, groups, unique_groups, Rs2, 
          eta2, rcut,  natoms, descr_size,two_body_decay)

    return rep
    
def generate_local_two_body_fchl_acsf(idx, nuclear_charges, coordinates, elements = [1,6,7,8,16],
        nRs2=24,  eta2=0.32, rcut=8.0, two_body_decay=1.8, 
        pad=False):

    Rs2 = np.linspace(0, rcut, 1+nRs2)[1:]

    n_elements = len(elements)
    natoms = len(coordinates)

    descr_size = n_elements * nRs2 
    rep = fgenerate_local_two_body_fchl_acsf(idx, coordinates, nuclear_charges, elements, Rs2, 
          eta2, rcut,  natoms, descr_size,two_body_decay)

    return rep
    
def generate_local_fchl_acsf(idx, nuclear_charges, coordinates, elements = [1,6,7,8,16],
        nRs2=24, nRs3=20, nFourier=1, eta2=0.32, eta3=2.7, zeta=np.pi, rcut=8.0, acut=8.0,
        two_body_decay=1.8, three_body_decay=0.57, three_body_weight=13.4,
        pad=False):

    Rs2 = np.linspace(0, rcut, 1+nRs2)[1:]
    Rs3 = np.linspace(0, acut, 1+nRs3)[1:]

    Ts = np.linspace(0, np.pi, 2*nFourier)
    n_elements = len(elements)
    natoms = len(coordinates)

    descr_size = n_elements * nRs2 + (n_elements * (n_elements + 1)) * nRs3* nFourier

    # Normalization constant for three-body 
    three_body_weight = np.sqrt(eta3/np.pi) * three_body_weight

    rep = fgenerate_local_fchl_acsf(idx, coordinates, nuclear_charges, elements, Rs2, Rs3,
            Ts, eta2, eta3, zeta, rcut, acut, natoms, descr_size,
            two_body_decay, three_body_decay, three_body_weight)

    return rep
