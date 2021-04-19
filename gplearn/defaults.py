#!/usr/bin/env python
# coding: utf-8
# Author Xinyu Li

"Some default settings used in testing adsorbate-slab complex stablity"

from ase.data import covalent_radii as cradii
from ase.data import atomic_numbers as an
from ase import Atoms, Atom
def get_median_energies(adsorbate):
    assert adsorbate in ["H", "C", "N", 'O', "S", 'H1', 'CO']
    if adsorbate == "H":
        median = {"Ag":0.409,
                  "Al":0.351,
                  "Au":0.384,
                    "Bi":0.979,
                    "Cd":0.855,
                    "Co":-0.177,
                    "Cr":-0.877,
                    "Cu":0.022,
                    "Fe":-0.644,
                    "Ga":0.666,
                    "Hf":-0.866,
                    'Hg':0.799,
                    "In":0.907,
                    "Ir":-0.140,
                    "La":-0.545,
                    "Mo":-0.775,
                    'Mn':-0.454,
                    "Nb":-0.873,
                    'Ni':-0.216,
                    "Os":-0.199,
                    "Pb":0.806,
                    "Pd":-0.208,
                    "Pt":-0.216,
                    "Re":-0.415,
                    "Rh":-0.213,
                    "Ru":-0.323,
                    "Sc":-0.953,
                    "Sn":0.780,
                    "Ta":-0.952,
                    "Tc":-0.454,
                    "Ti":-0.956,
                    "Tl":0.987,
                    "V":-0.954,
                    "W":-0.863,
                    "Y":-0.903,
                    "Zn":0.580,
                    "Zr":-0.889}
    elif adsorbate =='N':
        median = {"Ag":3.207,
                    "Al":-0.932,
                    "Au":3.030,
                    "Bi":1.654,
                    "Co":-0.065,
                    "Cr":-2.285,
                    "Cu":1.606,
                    "Fe":-1.497,
                    "Ga":1.029,
                    "Hf":-2.926,
                    "In":1.820,
                    "Ir":0.340,
                    "La":-1.978,
                    "Mo":-1.906,
                    "Nb":-2.629,
                    "Os":0.072,
                    "Pb":2.036,
                    "Pd":1.096,
                    "Pt":0.994,
                    "Re":-0.888,
                    "Rh":0.265,
                    "Ru":-0.092,
                    "Sc":-2.439,
                    "Sn":1.358,
                    "Ta":-2.771,
                    "Tc":-0.869,
                    "Ti":-2.604,
                    "Tl":2.925,
                    "V":-2.680,
                    "W":-2.321,
                    "Y":-2.285,
                    "Zn":2.261,
                    "Zr":-2.602,
                    'Cd':2.261,
                    'Hg':2.261,
                    'Mn':-0.879,
                    'Ni':1.045}
    elif adsorbate == 'O':
        median = {"Ag":2.051,
                    "Al":-1.589,
                    "Au":2.619,
                    "Bi":0.986,
                    "Cd":0.965,
                    "Co":-0.201,
                    "Cr":-2.499,
                    "Cu":0.927,
                    "Fe":-1.072,
                    "Ga":0.358,
                    "Hf":-3.325,
                    "Hg":1.846,
                    "In":0.727,
                    "Ir":0.907,
                    "La":-2.912,
                    "Mo":-1.927,
                    "Nb":-2.977,
                    "Os":0.287,
                    "Pb":0.903,
                    "Pd":1.499,
                    "Pt":1.804,
                    "Re":-0.962,
                    "Rh":0.598,
                    "Ru":0.046,
                    "Sc":-3.692,
                    "Sn":0.702,
                    "Ta":-3.040,
                    "Tc":-0.820,
                    "Ti":-3.254,
                    "Tl":1.399,
                    "V":-3.167,
                    "W":-2.012,
                    "Y":-3.628,
                    "Zn":0.193,
                    "Zr":-3.136,
                    'Mn':-0.879,
                    'Ni':1.045}
    elif adsorbate == 'C':
        median = {"Ag":5.667,
                    "Al":2.956,
                    "Au":4.897,
                    "Cd":5.411,
                    "Co":2.293,
                    "Cr":1.071,
                    "Cu":4.163,
                    "Fe":2.313,
                    "Ga":4.109,
                    "Hf":1.118,
                    "In":4.909,
                    "Ir":2.245,
                    "La":2.210,
                    "Mo":1.135,
                    "Nb":0.911,
                    "Os":2.302,
                    "Pb":5.196,
                    "Pd":2.559,
                    "Pt":2.306,
                    "Re":1.640,
                    "Rh":2.220,
                    "Ru":2.166,
                    "Sc":1.980,
                    "Sn":4.464,
                    "Ta":0.752,
                    "Tc":1.716,
                    "Ti":1.310,
                    "Tl":5.911,
                    "V":1.009,
                    "W":0.722,
                    "Y":2.209,
                    "Zn":5.302,
                    "Zr":1.283,
                    'Bi':1.654,
                    'Hg':2.261,
                    'Mn':-0.879,
                    'Ni':1.045}
    elif adsorbate == 'S':
        median = {"Ag":-0.243,
                    "Al":-1.492,
                    "Au":-0.029,
                    "Bi":-0.360,
                    "Cd":-0.354,
                    "Co":-1.869,
                    "Cr":-3.046,
                    "Cu":-1.035,
                    "Fe":-3.049,
                    "Ga":-0.485,
                    "Hf":-3.239,
                    "Hg":1.060,
                    "In":-0.232,
                    "Ir":-1.920,
                    "La":-2.948,
                    "Mo":-2.903,
                    "Nb":-3.332,
                    "Os":-1.726,
                    "Pb":-0.178,
                    "Pd":-1.276,
                    "Pt":-1.337,
                    "Re":-1.774,
                    "Rh":-1.908,
                    "Ru":-2.043,
                    "Sc":-3.448,
                    "Sn":-0.453,
                    "Ta":-3.340,
                    "Tc":-2.430,
                    "Ti":-3.174,
                    "Tl":0.044,
                    "V":-3.323,
                    "W":-3.048,
                    "Y":-3.312,
                    "Zn":-0.803,
                    "Zr":-3.137,
                    'Mn':-0.879,
                    'Ni':1.045}
    elif adsorbate == "H1":
        median = {'Hf': -0.9554726949999996,
                 'Ta': -0.5605048150000038,
                 'Mn': -0.40109448499998956,
                 'Si': -2.658273614999987,
                 'Co': -0.4352923150000003,
                 'Ir': -0.6032846749999998,
                 'Au': 0.41992627499999946,
                 'Sr': 0.49327162499999977,
                 'Zr': -0.5287337749999943,
                 'V': -0.24086023500001597,
                 'Te': 0.7705334449999959,
                 'In': 0.5671323350000015,
                 'Se': 0.8964256749999939,
                 'Cd': 0.6149637050000001,
                 'W': -0.5097040350000293,
                 'Ca': 0.4206683049999973,
                 'Ni': 0.0699429150000106,
                 'Rh': -0.28172470499999,
                 'Y': -0.34660050500000006,
                 'Tc': -0.49364522500002517,
                 'Sn': 1.6057355050000033,
                 'N': 2.3622441549999826,
                 'Nb': -0.568394285000013,
                 'Ti': -0.9902275549999922,
                 'Ga': 0.3565817949999981,
                 'Os': -0.5115438850000111,
                 'Sc': -0.5696183150000018,
                 'Bi': 0.452593264999988,
                 'Pt': -0.32401606499999724,
                 'Mo': -0.6887495549999891,
                 'Ag': 0.4098794449999992,
                 'Re': -0.7609179049999875,
                 'Fe': -0.4020598449999988,
                 'K': 0.274312895,
                 'Sb': 0.2503862249999975,
                 'P': -0.09510602500000909,
                 'Zn': 0.43893605499999966,
                 'Pb': 1.137056694999996,
                 'S': -0.09510602500000909,
                 'As': 1.371823705000009,
                 'Ru': -0.4294440850000023,
                 'Na': 0.34904644500000037,
                 'H': 0.0000000,
                 'Cr': -0.5752514849999693,
                 'Cu': 0.13348939500001,
                 'Ge': -1.665093335000011,
                 'Pd': -0.3190505449999983,
                 'Al': 0.17089637499999588,
                 'C': -0.09510602500000909}
    elif adsorbate == "CO":
        median = {'Hf': -1.0655779600000042,
                 'Ta': -1.0093796099999874,
                 'Mn': -1.327984919999997,
                 'Si': -0.35570180999999046,
                 'Co': -1.7767136200000007,
                 'Ir': -1.3049350699999902,
                 'Au': 0.09079323999999112,
                 'Sr': -0.03218879000000108,
                 'Zr': -1.0545259099999864,
                 'V': -0.01789497999998524,
                 'Te': 0.004451789999999178,
                 'In': 0.020582339999998922,
                 'Se': -0.5218858700000037,
                 'Cd': 0.009012569999997666,
                 'Ca': -0.29352835000000077,
                 'W': -1.6650268799999903,
                 'Ni': -1.7256401400000012,
                 'Rh': -1.6994628199999706,
                 'Y': -0.7142085199999979,
                 'Tc': -1.5684091000000446,
                 'Sn': 1.304411080000003,
                 'Cs': 0.001432219999999873,
                 'N': -0.014919089999997581,
                 'Nb': -1.2030138800000056,
                 'Ti': -1.8624212900000199,
                 'Ga': 0.1161423000000017,
                 'Os': -1.9894173599999956,
                 'Sc': -0.8044173599999933,
                 'Bi': -0.012824490000019395,
                 'Pt': -1.3821357100000196,
                 'Ag': 0.035632379999997355,
                 'Mo': -1.5550309699999811,
                 'Re': -1.5154386199999887,
                 'Fe': -1.7144563600000051,
                 'K': -0.059992230000000646,
                 'Sb': 0.03264844000001155,
                 'P': -0.4795591200000242,
                 'Zn': -0.09162822000000048,
                 'Pb': -0.010076760000002238,
                 'S': -0.4795591200000242,
                 'As': 0.01472578000000091,
                 'Ru': -1.6344762099999937,
                 'Na': -0.08089534999999692,
                 'H': -0.005110469999999978,
                 'Cr': -1.2233117199999608,
                 'Cu': -0.43935864999999374,
                 'Ge': -1.4305189799999933,
                 'Pd': -1.1948809700000051,
                 'Al': -0.14126516999998806}
    return median
    
elements_list = ['Ag', 'Al', 'Au', 'Bi', 'Cd', 'Co',
                 'Cr', 'Cu', 'Fe', 'Ga', 'Hf', 'Hg', 'In', 'Ir',
                 'La', 'Mo', 'Mn', 'Nb', 'Ni', 'Os', 'Pb', 'Pd',
                 'Pt', 'Re', 'Rh', 'Ru', 'Sc', 'Sn', 'Ta', 'Tc',
                 'Ti', 'Tl', 'V', 'W', 'Y', 'Zn', 'Zr']

elements_number = [an[n] for n in elements_list]

natcat_h_alloy_num = [22, 23, 24, 25, 26, 27, 28, 29, 41, 42, 44, 45, 46, 47, 74, 75,
       76, 77, 78, 79] # Only contain transition metallics
natcat_co_alloy_num = [22, 23, 24, 25, 26, 27, 28, 29, 41, 42, 44, 45, 46, 47, 74,
       75, 76, 77, 78, 79] # Only contain transition metallics
natcat_metallics_num = [13, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 41, 42, 44, 45,
       46, 47, 48, 49, 50, 74, 75, 76, 77, 78, 79, 82]

natcat_h_ele_num = [ 1,  6,  7, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                       28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                       48, 49, 50, 51, 52, 72, 73, 74, 75, 76, 77, 78, 79, 82, 83]
natcat_co_ele_num = [ 1,  7, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                   29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                   49, 50, 51, 52, 55, 72, 73, 74, 75, 76, 77, 78, 79, 82, 83]
natcat_co_ele_num = [ 1,  6, 7, 8, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                   29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                   49, 50, 51, 52, 55, 72, 73, 74, 75, 76, 77, 78, 79, 82, 83]

def defaults_ads_height(site, meancradii, adsorbate = 'H'):
    if adsorbate == 'H':
        if site == 'top':
            ads_height = 0.90 * meancradii
        elif site =='bridge':
            ads_height = 1.1
        elif site == 'hole':
            ads_height = 1.0
    elif adsorbate =='N':
        if site == 'top':
            ads_height = 0.80 * meancradii
        elif site =='bridge':
            ads_height = 1.0
        elif site == 'hole':
            ads_height = 1.0
    elif adsorbate =='NH':
        if site == 'top':
            ads_height = 0.80 * meancradii
        elif site =='bridge':
            ads_height = 1.20
        elif site == 'hole':
            ads_height = 1.10
    elif adsorbate =='O':
        if site == 'top':
            ads_height = 0.80 * meancradii
        elif site =='bridge':
            ads_height = 1.15
        elif site == 'hole':
            ads_height = 1.10
    elif adsorbate =='OH':
        if site == 'top':
            ads_height = 0.84* meancradii
        elif site =='bridge':
            ads_height = 1.45
        elif site == 'hole':
            ads_height = 1.35
    elif adsorbate =='C':
        if site == 'top':
            ads_height = 0.70 * meancradii
        elif site =='bridge':
            ads_height = 1.0
        elif site == 'hole':
            ads_height = 1.0
    elif adsorbate =='CH':
        if site == 'top':
            ads_height = 0.75 * meancradii
        elif site =='bridge':
            ads_height = 1.15
        elif site == 'hole':
            ads_height = 1.15
    elif adsorbate =='CH2':
        if site == 'top':
            ads_height = 0.77 * meancradii
        elif site =='bridge':
            ads_height = 1.45
        elif site == 'hole':
            ads_height = 1.45
    elif adsorbate =='CH3':
        if site == 'top':
            ads_height = 0.80 * meancradii
        elif site =='bridge':
            ads_height = 1.60
        elif site == 'hole':
            ads_height = 1.60
    elif adsorbate =='S':
        if site == 'top':
            ads_height = 0.80 * meancradii
        elif site =='bridge':
            ads_height = 1.70
        elif site == 'hole':
            ads_height = 1.70
    elif adsorbate =='SH':
        if site == 'top':
            ads_height = 0.90 * meancradii
        elif site =='bridge':
            ads_height = 1.90
        elif site == 'hole':
            ads_height = 1.90
    return ads_height

PTP = {
         1  :[1,1] ,2:  [1,8]#Row1

        ,3  :[2,1] ,4:  [2,2]#Row2\
        ,5  :[2,3] ,6:  [2,4] ,7  :[2,5] ,8  :[2,6] ,9  :[2,7] ,10 :[2,8]\

        ,11 :[3,1] ,12: [3,2]#Row3\
        ,13 :[3,3] ,14: [3,4] ,15 :[3,5] ,16 :[3,6] ,17 :[3,7] ,18 :[3,8]\

        ,19 :[4,1] ,20: [4,2]#Row4\
        ,31 :[4,3] ,32: [4,4] ,33 :[4,5] ,34 :[4,6] ,35 :[4,7] ,36 :[4,8]\
        ,21 :[4,9] ,22: [4,10],23 :[4,11],24 :[4,12],25 :[4,13],26 :[4,14],27 :[4,15],28 :[4,16],29 :[4,17],30 :[4,18]\

        ,37 :[5,1] ,38: [5,2]#Row5\
        ,49 :[5,3] ,50: [5,4] ,51 :[5,5] ,52 :[5,6] ,53 :[5,7] ,54 :[5,8]\
        ,39 :[5,9] ,40: [5,10],41 :[5,11],42 :[5,12],43 :[5,13],44 :[5,14],45 :[5,15],46 :[5,16],47 :[5,17],48 :[5,18]\

        ,55 :[6,1] ,56: [6,2]#Row6\
        ,81 :[6,3] ,82: [6,4] ,83 :[6,5] ,84 :[6,6] ,85 :[6,7] ,86 :[6,8]
               ,72: [6,10],73 :[6,11],74 :[6,12],75 :[6,13],76 :[6,14],77 :[6,15],78 :[6,16],79 :[6,17],80 :[6,18]\
        ,57 :[6,19],58: [6,20],59 :[6,21],60 :[6,22],61 :[6,23],62 :[6,24],63 :[6,25],64 :[6,26],65 :[6,27],66 :[6,28],67 :[6,29],68 :[6,30],69 :[6,31],70 :[6,32],71 :[6,33]\

        ,87 :[7,1] ,88: [7,2]#Row7\
        ,113:[7,3] ,114:[7,4] ,115:[7,5] ,116:[7,6] ,117:[7,7] ,118:[7,8]\
               ,104:[7,10],105:[7,11],106:[7,12],107:[7,13],108:[7,14],109:[7,15],110:[7,16],111:[7,17],112:[7,18]\
        ,89 :[7,19],90: [7,20],91 :[7,21],92 :[7,22],93 :[7,23],94 :[7,24],95 :[7,25],96 :[7,26],97 :[7,27],98 :[7,28],99 :[7,29],100:[7,30],101:[7,31],101:[7,32],102:[7,14],103:[7,33]}

def group(nuclear_charge):
    """
    the funtion used to get the group id
    intput:
    an   int
    elemenet number
    """
    return PTP[nuclear_charge][1]

def period(nuclear_charge):
    """
    the funtion used to get the period number
    intput:
    an   int
    elemenet number
    """
    return PTP[nuclear_charge][0]

def default_adsorbates():
    '''
    A dictionary whose keys are the simple string names of adsorbates and whos
    values are their corresponding `ase.Atoms` objects. When making new entries
    for this dictionary, we recommend "pointing" the adsorbate upwards in the
    z-direction.
    '''
    adsorbates = {}
    adsorbates[''] = Atoms()

    # Uranium is a place-holder for an adsorbate
    adsorbates['U'] = Atoms('U')

    # We put some of these adsorbates closer to the slab to help them adsorb
    # onto the surface
    adsorbates['H'] = Atoms('H')
    adsorbates['N'] = Atoms('N')
    adsorbates['O'] = Atoms('O')
    adsorbates['C'] = Atoms('C')
    adsorbates['S'] = Atoms('S')
    

    # For diatomics (and above), it's a good practice to manually relax the gases
    # and then see how far apart they are. Then put first atom at the origin, and
    # put the second atom directly above it.
    """
    adsorbates['CH'] = Atoms('CH', positions=[[0., 0., 0.],
                                              [0., 0., 1.1]])
    adsorbates['NH'] = Atoms('NH', positions=[[0., 0., 0.],
                                              [0., 0., 1.05]])
    adsorbates['OH'] = Atoms('OH', positions=[[0., 0., 0.],
                                              [0., 0., 0.95]])
    adsorbates['SH'] = Atoms('OH', positions=[[0., 0., 0.],
                                              [0., 0., 1.25]])
    
    adsorbates['CH2'] = Atoms('CH2', positions = [[ 0.        ,  0.        ,  0.        ],
                                          [ 0.62999996,  0.62999996,  0.62999996],
                                          [-0.62999996, -0.62999996,  0.62999996]])
    adsorbates['CH3'] = Atoms('CH3', positions =[[ 0.        ,  0.        ,  0.        ],
                                           [ 1.08999992,  0.        ,  0.62999996],
                                           [-0.44999997,  0.88949994,  0.62999996],
                                           [-0.44999997, -0.88949994,  0.62999996]])
    """
    adsorbates['CH'] = Atoms("C")
    adsorbates['CH2'] = Atoms("C")
    adsorbates['CH3'] = Atoms("C")
    adsorbates['NH'] = Atoms("N")
    adsorbates['OH'] = Atoms("O")
    adsorbates['SH'] = Atoms("S")
    return adsorbates

        
