#!/usr/bin/env python
# coding: utf-8
# GNU Lesser General Public License v3.0

"""
This source code is licensed under the GNU Lesser General Public License v3.0
found in the LICENSE file in the root directory of this source tree.

---
This code borrows heavily from
https://github.com/SUNCAT-Center/CatHub/blob/master/cathub/classification.py

----

Modifications made to enable site rebuilt
"""

import sys
import os
import numpy as np
import ase
from ase.geometry import get_distances, wrap_positions
from ase.db import connect
from ase.io import write
from ase.build import add_adsorbate
from scipy.spatial import Voronoi
from collections import OrderedDict
from ase.visualize import view
from copy import copy, deepcopy
from ase import Atoms, Atom
from ase.data import covalent_radii as cradii
from ase.data import atomic_numbers as an
from .defaults import defaults_ads_height, default_adsorbates

class SiteClassification:
    """ 
    Determine surface reconstruction (True/False) and adsorption
    site for chemisorption.
    A: ASE Atoms object
        Initial surface structure. Adsorbate atoms must be last indices
    B: ASE Atoms object
        Final surface structure. Adsorbate atoms must be last indices
    natoms_top_layer: int
        Number of atoms in top layer of slab
    natoms_slab: int
        Number of atoms in slab
    """

    def __init__(self, B, natoms_top_layer=4, natoms_slab=12, A=None):
        self.ntop = natoms_top_layer
        self.nslab = natoms_slab
        self.B = B
        self.initialA = A.copy()
        self.initialB = B.copy()
        # Check if multiatomic adsorbates dissociates on surface
        self.dissociated = self.check_dissociated()

        # Only keep the adsorbate closest to surface
        B = self.remove_extra_atoms(B)

        # Sort layers according to z-position
        layer_indices = np.argsort(B[:-1].positions[:, 2])
        self.B = B[:-1][layer_indices] + B[-1]
        
        if A:
            A = self.remove_extra_atoms(A)
            self.A = A
            self.A = A[:-1][layer_indices] + A[-1]

    def get_info(self):
        """Return surface reconstruction as well as primary and
        secondary adsorption site labels"""

        reconstructed = self.is_reconstructed()

        site, site_type = self.get_site()

        return reconstructed, site, site_type

    def check_dissociated(self, cutoff=1.2):
        """Check if adsorbate dissociates"""
        dissociated = False
        if not len(self.B) > self.nslab + 1:  # only one adsorbate
            return dissociated

        adsatoms = [atom for atom in self.B[self.nslab:]]
        ads0, ads1 = set(atom.symbol for atom in adsatoms)
        bond_dist = get_ads_dist(self.B, ads0, ads1)

        Cradii = [cradii[atom.number]
                  for atom in [ase.Atom(ads0), ase.Atom(ads1)]]
        bond_dist0 = sum(Cradii)

        if bond_dist > cutoff * bond_dist0:
            print('DISSOCIATED: {} Ang > 1.2 * {} Ang'
                  .format(bond_dist, bond_dist0))
            dissociated = True

            return dissociated

    def remove_extra_atoms(self, slab):
        if not len(slab) > self.nslab + 1:
            return slab

        adsatoms = [atom for atom in slab[self.nslab:]]
        adsindex = np.argmin([atom.position[2]
                              for atom in adsatoms])

        del slab[[atom.index + self.nslab for atom in adsatoms
                  if not atom.index == adsindex]]

        return slab

    def is_desorbed(self):
        desorbed = False
        D, D_len = get_distances(self.B.positions, cell=self.B.cell, pbc=True)

        indexM = np.argmin(D_len[-1, :-1])
        dist_S = D_len[-1, indexM]

        if dist_S > (cradii[self.B[-1].number] +
                     cradii[self.B[indexM].number]) * 2:
            print('DESORBED FROM SLAB')
            desorbed = True
        return desorbed

    def is_reconstructed(self, xy_cutoff=0.3, z_cutoff=0.4):
        """Compare initial and final slab configuration 
        to determine if slab reconstructs during relaxation
        xy_cutoff: Allowed xy-movement is determined from 
                   the covalent radii as:
                   xy_cutoff * np.mean(cradii)
        z_cutoff:  Allowed z-movement is determined as
                   z_cutoff * cradii_i
        """

        assert self.A, \
            'Initial slab geometry needed to classify reconstruction'

        # remove adsorbate
        A = self.A[:-1].copy()
        B = self.B[:-1].copy()

        # Order wrt x-positions
        x_indices = np.argsort(A.positions[:, 0])
        A = A[x_indices]
        B = B[x_indices]
        a = A.positions
        b = B.positions

        allowed_z_movement = z_cutoff * cradii[A.get_atomic_numbers()]
        allowed_xy_movement = \
            xy_cutoff * np.mean(cradii[A.get_atomic_numbers()])

        D, D_len = get_distances(p1=a, p2=b, cell=A.cell, pbc=True)
        d_xy = np.linalg.norm(np.diagonal(D)[:2], axis=0)
        d_z = np.diagonal(D)[2:][0]

        cond1 = np.all(d_xy < allowed_xy_movement)
        cond2 = np.all([d_z[i] < allowed_z_movement[i]
                        for i in range(len(a))])

        if cond1 and cond2:  # not reconstructed
            return False
        else:
            return True

    def is_subsurface(self):
        pos0 = self.B.positions[:-1][:, 2]
        pos1 = self.B.positions[-1][2]
        metal_covalent_radii = cradii[self.B.get_atomic_numbers()[:-1]]
        ads_covalent_radii = cradii[self.B.get_atomic_numbers()[-1]]

        if np.any([(pos0[i] - pos1) > 0.5 * metal_covalent_radii[i]
                   for i in range(len(pos0))]):
            return True
        else:
            return False

    def get_site_dict(self, ads_pos):
        """Get dictionary with high symmetry sites close to adsorbate
        position. 
        Top sites: Optained from the atomic positions of
                   the top layer of the slab. 
        Bridge sites: The average position of atomic pairs
        Hollow sites: Optained as the Voronoi vertices
        4-fold sites: Assigned when Voronoi vertives overlap
        """

        # Get top layer
        C = self.initialB[-self.ntop - 1:-1]
        IC = self.initialA[-self.ntop - 1:-1]
        
        C = get_real_space_positions(IC, C)
        
        SC = C * (3, 3, 1)
        ISC = IC * (3, 3, 1)
        self.__SC = SC
        self.__ISC = ISC
        
        D = OrderedDict()
        cell = C.get_cell()
        top_dist = 1
        """Top sites: atomic positions """
        for i, atom in \
            enumerate([atom for atom in SC if
                       np.linalg.norm(atom.position[:2] - ads_pos)
                       < top_dist]):
            k = 'top_site-{}'.format(i)
            D[k] = {}
            D[k]['pos'] = atom.position
            D[k]['sym'] = atom.symbol
            D[k]['ipos'] = ISC[atom.index].position

        """Bridge sites: bewteen atomic pairs """
        bridge_dist = 0.5 * cell[0][0]
        b_index = 0
        for atomi in [atom for atom in SC if
                      np.linalg.norm(atom.position[:2] - ads_pos)
                      < bridge_dist]:
            pos1 = atomi.position
            # Atom bairs should be close to each other and close to
            # adsorbate
            for atom in \
                [atom for atom in SC if np.linalg.norm(atom.position - pos1)
                 < 1.5 * bridge_dist and
                 np.linalg.norm(atom.position[:2] - ads_pos) < bridge_dist
                 and not (atom.position == pos1).all()]:
                pos2 = atom.position
                # bridge site is average position
                bridge_pos = 0.5 * (pos1 + pos2)
                k = 'bridge_site-{}'.format(b_index)
                D[k] = {}
                D[k]['pos'] = bridge_pos
                D[k]['sym'] = atomi.symbol + '_' + atom.symbol
                D[k]['ipos'] = 0.5 * (ISC[atomi.index].position + ISC[atom.index].position)
                b_index += 1 

        """Hollow sites: Voronoi vertices """
        hollow_dist = 1
        vor = Voronoi(SC.positions[:, :2])
        vertices = vor.vertices
        # map vertices close to adsorbate
        close_v = [v for v in vertices if np.linalg.norm(v - ads_pos[:2])
                   < hollow_dist]
        h_index = 1
        ff_index = 1
        for v in close_v:
            # Check if vertices overlap
            close_v_v = [v0 for v0 in close_v if np.linalg.norm(v0 - v) < 0.5]
            if len(close_v_v) > 1:
                v_mean = np.mean(close_v_v, axis=0)
                k = '4fold_{}'.format(ff_index)
                D[k] = {}
                D[k]['pos'] = v_mean

                # Delete bridge sites overlapping with 4fold
                for key in [key for key in list(D.keys()) if 'bridge' in key]:
                    bridge_pos = D[key]['pos']
                    if np.linalg.norm(bridge_pos[:2] - v) < 0.3:
                        del D[key]

                ffold_max_dist = sorted([np.linalg.norm(v_mean - m[:2])
                                         for m in SC.positions])[4]
                symb = [atom.symbol for atom in SC if
                        np.linalg.norm(v_mean - atom.position[:2])
                        < ffold_max_dist]

                D[k]['sym'] = '_'.join(symb)
                D[k]['ipos'] = 1/4 * np.array([ISC[idx].position for idx, atom in enumerate(SC) if
                        np.linalg.norm(v_mean - atom.position[:2])
                        < ffold_max_dist]).sum(axis = 0)
                ff_index += 1
            else:  # Regular hollow site
                k = 'hollow_site-{}'.format(h_index)
                D[k] = {}
                D[k]['pos'] = v
                hollow_max_dist = sorted([np.linalg.norm(v - m[:2])
                                          for m in SC.positions])[3]
                symb = [atom.symbol for atom in SC
                        if np.linalg.norm(v - atom.position[:2])
                        < hollow_max_dist]
                D[k]['sym'] = '_'.join(symb)
                
                D[k]['ipos'] = 1/3 * np.array([ISC[idx].position for idx, atom in enumerate(SC) if
                        np.linalg.norm(v - atom.position[:2])
                        < hollow_max_dist]).sum(axis = 0)
                D[k]['idx'] = [idx for idx, atom in enumerate(SC) if
                        np.linalg.norm(v - atom.position[:2])
                        < hollow_max_dist]
                h_index += 1

        return D

    def get_subsurface_layer(self):
        return self.B[np.argsort(self.B.positions[:, 2])][-self.ntop * 2 - 1:
                                                          -self.ntop - 1]

    def get_under_bridge(self):
        """Return element closest to the adsorbate in the subsurface layer"""
        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.get_subsurface_layer() * (3, 3, 1)
        dis = self.B.cell[0][0] * 2

        ret = None

        for ele in C:
            new_dis = np.linalg.norm(ads_pos - ele.position)
            if new_dis < dis:
                dis = new_dis
                ret = ele.symbol

        return ret

    def get_under_hollow(self):
        """ Return HCP if an atom is present below the adsorbate in the 
        subsurface layer and FCC if not"""
        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.get_subsurface_layer() * (3, 3, 1)

        ret = 'FCC'
        if np.any([np.linalg.norm(ads_pos[:2] - ele.position[:2]) < 0.5 *
                   cradii[ele.number] for ele in C]):
            ret = 'HCP'

        return ret

    def get_site(self, plot_voronoi_sites=False):
        """Return primaty (top, bridge, hollow, 4fold) and
        secondary (chemical elements in close environment) site designation"""

        if self.dissociated:
            return 'dissociated', ''

        if self.is_desorbed():
            return 'desorbed', ''

        if self.is_subsurface():
            return 'subsurface', ''

        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.B.copy() * (3, 3, 1)

        # Use top layer and adsorbate to map sites
        Dict = self.get_site_dict(ads_pos[:2])
        primary_site = None
        dis = self.B.get_cell()[0][0]
        Kind = None

        values = [np.linalg.norm(ads_pos[:2] - d['pos'][:2])
                  for d in list(Dict.values())]
        if len(values) == 0:
            return 'N/A', ''
        idx = np.argmin(values)
        
        dis = values[idx]
        kind = list(Dict.keys())[idx]
        self.ipos = Dict[kind]['ipos']
        primary_site = kind.split('_')[0]

        if plot_voronoi_sites:  # View sampled sites
            X = self.B.copy()
            X = X * (3, 3, 1)
            del X[-1]
            for pos in Dict.values():
                add_adsorbate(X, 'X', position=(pos['pos'][:2]), height=0.2)
            view(X)

        if primary_site == 'top':
            site_type = Dict[kind]['sym']

        if primary_site == 'bridge':
            site_type = Dict[kind]['sym'] + '|' + self.get_under_bridge()

        elif primary_site == 'hollow':
            site_type = Dict[kind]['sym'] + '|' + self.get_under_hollow()

        elif primary_site == '4fold':
            site_type = Dict[kind]['sym']

        if dis > 0.5:
            primary_site += '-tilt'
            print('Warning: A strong site match could not be found!')
            print('  structure labeled as {}'.format(primary_site))

        return primary_site, site_type
    
    def rebuild_initial(self, adsorbate = 'H'):
        primary_site, site_type = self.get_site()
        buttom = [atom for atom in self.A[:self.ntop]]
        Cradii_slab = np.min([cradii[atom.number] for atom in buttom])
        ads_height = cradii[Atom(adsorbate).number] + Cradii_slab
        ad_site = self.ipos
        
        if 'top' in primary_site:
            ad_site[-1] += 0.85 * ads_height
        elif 'bridge' in primary_site:
            ad_site[-1] += 0.75 * ads_height
        elif 'hollow' in primary_site:
            ad_site[-1] += 0.65 * ads_height
        adsorbate = Atoms(adsorbate, positions=[[0., 0., 0.]]) 
        adsorbate.translate(ad_site)
        rebuild = self.initialA[:self.nslab] + adsorbate
        rebuild.wrap()
        return rebuild
        
def get_real_space_positions(initial, final):
    initial_, final_ = initial.copy(), final.copy()
    assert len(initial) == len(final)
    natoms = len(initial)
    for idx in range(natoms):
        at = final_[[idx]]
        iat = initial_[[idx]]
        superat = at * (3, 3, 1)
        xycell = at.cell.array[[0, 1]]
        superat.positions -= xycell.sum(axis= 0)
        dist = np.linalg.norm(superat.positions - iat.positions, axis =1)
        if np.argmin(dist) != 4:
            final_.positions[[idx]] = superat.positions[np.argmin(dist)]
    return final_        
        


def get_ads_dist(atoms, ads0, ads1='H'):
    index0 = [i for i in range(len(atoms)) if atoms[i].symbol == ads0]
    index1 = [i for i in range(len(atoms)) if atoms[i].symbol == ads1]
    dist = []
    D, D_len = get_distances(atoms.positions[index0],
                             atoms.positions[index1],
                             atoms.cell, pbc=True)

    return np.max(D_len)

class SiteRebuild:
    """ 
    Reconstruct the initial configurations used in ScientificData 
    www.nature.com/articles/s41597-019-0080-z
    atoms: ASE Atoms object
        Initial surface structure. Adsorbate atoms must be last indices
    slab_type: str
        A1, L10, L12 
    metalA: str
        element name of metal A
    metalB: str
        element name of metal B
    adsorbate: str
        the name of adsorbate
    """
    

    def __init__(self, atoms, slab_type, metalA ='Ni', metalB = 'Y',  adsorbate="H"):
        assert slab_type in ['A1', 'L12', 'L10']
        self.atoms = atoms
        self.slab_type = slab_type
        self.A = metalA
        self.B = metalB
        self.adsorbate = adsorbate
        
        site_dict = {"A1": ["A", "A_A_A|HCP", "A_A|A", "A_A_A|FCC"],
                     "L10":["A_B|B", "B_B|A", "A_A|B", "A_B_B|HCP", "A", "B", "A_A_B|FCC", \
                            "A_B_B|FCC", "A_A_B|HCP", "A_B|A"],
                     "L12":["A_A_A|FCC", "A_A|A","A_B|A","A_A|B",
                               "A","B","A_A_B|FCC","A_A_A|HCP","A_A_B|HCP"]}
        self.slab = atoms[:12].copy()
        self.top_layer = self.slab[-4:].copy()
        self.pos = copy(self.slab.positions)
        self.get_slab_ads_height()
        self.site_dict = {}
        if self.slab_type == 'A1':
            self.get_A1_site_dict()      
        elif self.slab_type == 'L10': 
            self.get_L10_site_dict()
        elif self.slab_type == 'L12':
            self.get_L12_site_dict()
        
    def get_slab_ads_height(self):
        self.cradiiA = cradii[an[self.A]]
        if self.B is not None:
            self.cradiiB = cradii[an[self.B]]
        self.cradiiAds = cradii[an[self.adsorbate[0]]] # TODO: the trend should be different with adsorbate
        self.slab_height = np.max(self.pos[:, -1])

    def get_A1_site_dict(self):
        x_axis = copy(self.slab.cell.array[0, 0])
        sub_pos= copy(self.pos[4])
        buttom_pos = copy(self.pos[0]) 
        
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict['A'] = [0., 0., self.slab_height + defaults_ads_height('top', meancradii, adsorbate = self.adsorbate)]
        self.site_dict['A_A|A'] =  [1/4 * x_axis, 0, self.slab_height + defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)]
        self.site_dict['A_A_A|HCP'] = sub_pos
        self.site_dict['A_A_A|HCP'][2] =  self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        self.site_dict['A_A_A|FCC'] = buttom_pos
        self.site_dict['A_A_A|FCC'][2] =  self.slab_height +defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
    def get_L10_site_dict(self):
        pos = [atom.position for atom in self.top_layer if atom.symbol == self.A]
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict['A'] = [pos[0][0], pos[0][1], self.slab_height + defaults_ads_height('top', meancradii, adsorbate = self.adsorbate)]
        pos = [atom.position for atom in self.top_layer if atom.symbol == self.B]
        meancradii = self.cradiiB + self.cradiiAds
        self.site_dict['B'] = [pos[0][0], pos[0][1], self.slab_height + defaults_ads_height('top', meancradii, adsorbate = self.adsorbate)]
        
        # Get sites of HCP and FCC sites
        subslab = self.slab[-8:-4].copy()
        top = self.top_layer.copy()  
        stop = top*(3, 3, 1)
        xycell = top.cell.array[[0, 1]]
        stop.positions -= xycell.sum(axis= 0)
        
        Apos = np.array([atom.position for atom in subslab if atom.symbol == self.B])[0, :-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3] 
        self.site_dict['A_A_B|HCP'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = (2 * self.cradiiA + self.cradiiB + 3*self.cradiiAds)/3
        self.site_dict['A_A_B|HCP'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        Apos = np.array([atom.position for atom in subslab if atom.symbol == self.A])[0,:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_B_B|HCP'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = (self.cradiiA + 2 * self.cradiiB + 3*self.cradiiAds)/3
        self.site_dict['A_B_B|HCP'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        buttomslab = self.slab[:4].copy()
        Apos = np.array([atom.position for atom in buttomslab if atom.symbol == self.B])[0,:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_A_B|FCC'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = (2 * self.cradiiA + self.cradiiB + 3*self.cradiiAds)/3
        self.site_dict['A_A_B|FCC'][-1]=self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        Apos = np.array([atom.position for atom in buttomslab if atom.symbol == self.A])[0,:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_B_B|FCC'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = (self.cradiiA + 2 * self.cradiiB + 3*self.cradiiAds)/3
        self.site_dict['A_B_B|FCC'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        # Get Bridge site      
        top = self.top_layer.copy()  
        sub = self.slab[-8:-4].copy()
        stop = top*(3, 3, 1)
        ssub = sub*(3, 3, 1)
        A, B = self.A, self.B
        pos = top.positions
        xycell = top.cell.array[[0, 1]]
        stop.positions -= xycell.sum(axis= 0)
        ssub.positions -= xycell.sum(axis= 0)
        Apos = [atom.position for atom in top if atom.symbol == A][0]
        dist = np.linalg.norm(stop.positions - Apos, axis =1)
        sortidx = np.argsort(dist)[1:7]
        elements = np.array([atom.symbol for atom in stop])
        clostdist, closeele = dist[sortidx], elements[sortidx]
        for idx, ele in enumerate(closeele):
            if ele == A:
                if "A_A|B" not in self.site_dict.keys():
                    self.site_dict["A_A|B"] = 1/2 * (stop.positions[sortidx[idx]] + Apos)
                    meancradii = self.cradiiA + self.cradiiAds
                    self.site_dict['A_A|B'][2] =  self.slab_height + defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)
                    
            if ele == B:
                xy_pos = 1/2 * (stop[sortidx[idx]].position[:-1] + Apos[:-1])
                subdist = np.linalg.norm(ssub.positions[:, :-1] - xy_pos, axis =1)
                if ssub[np.argmin(subdist)].symbol == A:
                    if "A_B|A" not in self.site_dict.keys():
                        self.site_dict["A_B|A"] = np.zeros((3))
                        self.site_dict["A_B|A"][:-1] = xy_pos
                elif ssub[np.argmin(subdist)].symbol == B:
                    if "A_B|B" not in self.site_dict.keys():
                        self.site_dict["A_B|B"] = np.zeros((3))
                        self.site_dict["A_B|B"][:-1] = xy_pos
        meancradii = (self.cradiiA+self.cradiiB+2*self.cradiiAds)/2
        self.site_dict["A_B|A"][-1] = self.slab_height + defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)
        self.site_dict["A_B|B"][-1] = self.site_dict["A_B|A"][-1]
        Bpos = [atom.position for atom in top if atom.symbol == B][0]
        dist = np.linalg.norm(stop.positions - Bpos, axis =1)
        sortidx = np.argsort(dist)[1:7]
        elements = np.array([atom.symbol for atom in stop])
        clostdist, closeele = dist[sortidx], elements[sortidx]
        for idx, ele in enumerate(closeele):
            if ele == B:
                if "B_B|A" not in self.site_dict.keys():
                    self.site_dict["B_B|A"] = 1/2 * (stop.positions[sortidx[idx]] + Bpos)
                    meancradii = self.cradiiB + self.cradiiAds
                    self.site_dict['B_B|A'][-1] =  self.slab_height + defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)
                    break

    def get_L12_site_dict(self):
        
        # Get TOP site
        pos = [atom.position for atom in self.top_layer if atom.symbol == self.A]
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict['A'] = [pos[0][0], pos[0][1], 
                               self.slab_height + defaults_ads_height('top', meancradii, adsorbate = self.adsorbate)]
        pos = [atom.position for atom in self.top_layer if atom.symbol == self.B]
        meancradii = self.cradiiB + self.cradiiAds
        self.site_dict['B'] = [pos[0][0], pos[0][1], 
                               self.slab_height + defaults_ads_height('top', meancradii, adsorbate = self.adsorbate)]
        
        # Get FCC and HCP 
        subslab = self.slab[-8:-4].copy()
        top = self.top_layer.copy()  
        stop = top*(3, 3, 1)
        xycell = top.cell.array[[0, 1]]
        stop.positions -= xycell.sum(axis= 0)
        subslab = self.slab[-8:-4].copy()
        
        Apos = [atom.position for atom in subslab if atom.symbol == self.B][0][:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_A_A|HCP'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict['A_A_A|HCP'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        Apos = [atom.position for atom in subslab if atom.symbol == self.A][0][:-1]
        
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_A_B|HCP'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = (2 * self.cradiiA + self.cradiiB + 3*self.cradiiAds)/3
        self.site_dict['A_A_B|HCP'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        buttomslab = self.slab[:4].copy()
        Apos = [atom.position for atom in buttomslab if atom.symbol == self.A][0][:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_A_B|FCC'] = stop.positions[sortidx].mean(axis = 0)
        self.site_dict['A_A_B|FCC'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
        
        Apos = [atom.position for atom in buttomslab if atom.symbol == self.B][0][:-1]
        dist = np.linalg.norm(stop.positions[:, :-1] - Apos, axis =1)
        sortidx = np.argsort(dist)[:3]
        self.site_dict['A_A_A|FCC'] = stop.positions[sortidx].mean(axis = 0)
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict['A_A_A|FCC'][-1] = self.slab_height + defaults_ads_height('hole', meancradii, adsorbate = self.adsorbate)
                
        # Get Bridge site
        top = self.top_layer.copy()
        sub = self.slab[-8:-4].copy()
        stop = top*(3, 3, 1)
        ssub = sub*(3, 3, 1)
        A, B = self.A, self.B
        pos = top.positions
        xycell = top.cell.array[[0, 1]]
        stop.positions -= xycell.sum(axis= 0)
        ssub.positions -= xycell.sum(axis= 0)
        Apos = [atom.position for atom in top if atom.symbol == A][1]
        dist = np.linalg.norm(stop.positions - Apos, axis =1)
        sortidx = np.argsort(dist)[1:7]
        elements = np.array([atom.symbol for atom in stop])
        clostdist, closeele = dist[sortidx], elements[sortidx]
        for idx, ele in enumerate(closeele):
            if ele == B:
                if "A_B|A" not in self.site_dict.keys():
                    self.site_dict["A_B|A"] = 1/2 * (stop.positions[sortidx[idx]] + Apos)
            if ele == A:
                xy_pos = 1/2 * (stop[sortidx[idx]].position[:-1] + Apos[:-1])
                subdist = np.linalg.norm(ssub.positions[:, :-1] - xy_pos, axis =1)
                if ssub[np.argmin(subdist)].symbol == A:
                    if "A_A|A" not in self.site_dict.keys():
                        self.site_dict["A_A|A"] = np.zeros((3))
                        self.site_dict["A_A|A"][:-1] = xy_pos
                elif ssub[np.argmin(subdist)].symbol == B:
                    if "A_A|B" not in self.site_dict.keys():
                        self.site_dict["A_A|B"] = np.zeros((3))
                        self.site_dict["A_A|B"][:-1] = xy_pos
        meancradii = (self.cradiiA +self.cradiiB+2*self.cradiiAds)/2
        self.site_dict['A_B|A'][-1] = self.slab_height+ defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)
        meancradii = self.cradiiA + self.cradiiAds
        self.site_dict["A_A|B"][-1] = self.slab_height+ defaults_ads_height('bridge', meancradii, adsorbate = self.adsorbate)
        self.site_dict["A_A|A"][-1] =  self.site_dict["A_A|B"][-1]  
        
    def rebuild_initial(self, site):
        adsorbate = default_adsorbates()[self.adsorbate]
        adsite = self.site_dict[site]
        adsorbate.set_tags([1] * len(adsorbate))
        adsorbate.translate(adsite)
        rebuild = self.slab + adsorbate
        rebuild.wrap()
        return rebuild
    
    
def SiteLabel2AB(site_type, A, B=None):
    site_type = site_type.replace(A, "A")
    if B != None:
        site_type = site_type.replace(B, "B")
    split = site_type.split("|")
    if len(split) == 1:
        if len(split[0]) == 1:
            return site_type
        else:
            split1 = split[0].split("_")
            split1 = np.sort(split1)
            return "_".join(split1)
    else:
        split1 = split[0].split("_")
        split1 = np.sort(split1)
        return "_".join(split1) + "|" + split[1]
