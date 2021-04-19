import os
import sys
sys.path.append(os.path.dirname(os.getcwd()))


def test_L12():
    doc = {'slab_name': 'Ga3Nb',
         'reconstructed': 0,
         'site': 'hollow',
         'metalA': 'Ga',
         'metalB': 'Nb',
         'site_type': 'A_A_B|HCP',
         'fw_id': 25325,
         'energy': -0.7984391464923988,
         'adsorbate': 'S',
         'initial_configuration': {'atoms': {'atoms': [{'symbol': 'Ga',
             'position': [0.0, 0.0, 9.99999929233104],
             'tag': 0,
             'index': 0,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [-1.42958099883298, 2.47610682477359, 9.99999929233104],
             'tag': 0,
             'index': 1,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [1.42958099883298, 2.47610682477359, 9.99999929233104],
             'tag': 0,
             'index': 2,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Nb',
             'position': [2.85916189766597, 0.0, 9.99999929233104],
             'tag': 0,
             'index': 3,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [2.85916189766597, 1.6507378831824, 12.334495227126],
             'tag': 0,
             'index': 4,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [1.42958099883298, 4.12684480795598, 12.334495227126],
             'tag': 0,
             'index': 5,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [-1.42958099883298, 4.12684480795598, 12.334495227126],
             'tag': 0,
             'index': 6,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Nb',
             'position': [0.0, 1.6507378831824, 12.334495227126],
             'tag': 0,
             'index': 7,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [0.0, 3.30147586636478, 14.668991161921],
             'tag': 0,
             'index': 8,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [1.42958099883298, 0.825368941591198, 14.668991161921],
             'tag': 0,
             'index': 9,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Ga',
             'position': [4.28874289649896, 0.825368941591198, 14.668991161921],
             'tag': 0,
             'index': 10,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'Nb',
             'position': [2.85916189766597, 3.30147586636478, 14.668991161921],
             'tag': 0,
             'index': 11,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0},
            {'symbol': 'S',
             'position': [-0.714790449416495, 4.53952927875158, 17.1489909864191],
             'tag': 1,
             'index': 12,
             'charge': 0.0,
             'momentum': [0.0, 0.0, 0.0],
             'magmom': 0.0}],
           'cell': {'array': [[5.71832359533196, 0.0, 0.0],
             [-2.85916179766598, 4.95221364954719, 0.0],
             [0.0, 0.0, 24.668990254252]],
            'pbc': [True, True, True],
            '__ase_objtype__': 'cell'},
           'pbc': [True, True, True],
           'info': {},
           'constraints': [{'name': 'FixAtoms',
             'kwargs': {'indices': [0, 1, 2, 3, 4, 5, 6, 7]}}],
           'natoms': 13,
           'mass': 938.28611,
           'spacegroup': 'P1 (1)',
           'chemical_symbols': ['Ga', 'Nb', 'S'],
           'symbol_counts': {'Ga': 9, 'Nb': 3, 'S': 1},
           'volume': 698.5853508362653},
          'calc': {'calculator': {'module': 'ase.calculators.singlepoint',
            'class': 'SinglePointCalculator'}},
          'results': {'energy': -56182.8098533397,
           'forces': [[-0.0228221401008381, -0.045257843396535, -0.882248081653326],
            [-0.0505095787137775, 0.002901232807195, -0.882303617481609],
            [0.0502249575938232, 0.0290974315865886, -0.903591837436425],
            [-0.0477939795548274, -0.0275928220166041, 2.70814630221654],
            [0.0411975572842323, -0.0104937004424722, -0.0301271584026836],
            [0.0115902759407607, 0.040975156860781, -0.0299705782201613],
            [-0.00482698907499814, -0.00274928061036301, 0.0223896805491098],
            [0.133647484985538, 0.0771606056424026, 0.0842391097557388],
            [0.143595854447501, 0.392096831640525, 1.59808536784716],
            [0.411404787940517, -0.0717072958373152, 1.59811262154067],
            [-0.27973705238917, -0.161596661900296, 1.01310900503266],
            [0.288614814494774, 0.166620340367132, -2.03784248167813],
            [-0.674585992853535, -0.389453480480407, -2.25799858917987]],
           'fmax': 2.25799858917987},
          'reconstructed': 0,
          'fw_id': 25325,
          'site': 'bridge',
          'site_type': 'A_A|A',
          'coordination': 'Ga-Ga',
          'neighborcoord': ['Ga:Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb',
           'Ga:Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb'],
          'nextnearest': 'Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb-Nb',
          'nextnearestcoordination': 'Ga-Ga-Ga-Nb'},
         'coordination': 'Ga-Ga',
         'neighborcoord': ['Ga:Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb',
          'Ga:Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb'],
         'nextnearest': 'Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Ga-Nb-Nb-Nb-Nb',
         'nextnearestcoordination': 'Ga-Ga-Ga-Ga-Nb'}
    
    from gplearn.site_rebuilt import SiteRebuild, SiteClassification
    from gplearn.read_datasets import make_atoms_from_doc
    iatoms = make_atoms_from_doc(doc['initial_configuration'])
    sc = SiteRebuild(iatoms, 
                 'L12',
                 metalA =doc['metalA'], 
                 metalB =doc['metalB'],  
                 adsorbate='O')
    for site in ["A_A_A|FCC", "A_A|A","A_B|A","A_A|B",
                "A","B","A_A_B|FCC","A_A_A|HCP","A_A_B|HCP"]:
        rebuild = sc.rebuild_initial(site)
        scla = SiteClassification(rebuild, A = rebuild)
        got_site = scla.get_site()[1]
        got_site = got_site.replace(doc['metalA'], "A")
        got_site = got_site.replace(doc['metalB'], "B")
        if "|" not in site:
            assert site == got_site
        else:
            _A , _B = got_site.split("|")
            __A , __B = site.split("|")

            _A = sorted(_A.split("_"))
            _B = sorted(_B.split("_"))
            __A = sorted(__A.split("_"))
            __B = sorted(__B.split("_"))
            print(_A, __A)
            assert _A == __A
            assert _B == __B
    assert _B == __B
    
    
def test_L10():
    """
    """
    doc = {
         'slab_name': 'AuTa',
         'reconstructed': 0,
         'site': 'hollow',
         'metalA': 'Au',
         'metalB': 'Ta',
         'site_type': 'A_B_B|HCP',
         'fw_id': 50514,
         'energy': -2.5764415807934142,
         'adsorbate': 'S',
         'initial_configuration': {'atoms': {'atoms': [{'symbol': 'Ta',
         'position': [4.33894779294611, 0.878074837861365, 9.99999929233104],
         'tag': 0,
         'index': 0,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [1.40116670084377, 0.878074837861365, 9.99999929233104],
         'tag': 0,
         'index': 1,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [2.84748269849247, 3.40910025874854, 9.99999929233104],
         'tag': 0,
         'index': 2,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Ta',
         'position': [5.78526369059482, 3.40910025874854, 9.99999929233104],
         'tag': 0,
         'index': 3,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Ta',
         'position': [2.89263189529741, 1.70455007937427, 12.3922844230364],
         'tag': 0,
         'index': 4,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [5.83041288739975, 1.70455007937427, 12.3922844230364],
         'tag': 0,
         'index': 5,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Ta',
         'position': [4.33894779294611, 4.23557550026145, 12.3922844230364],
         'tag': 0,
         'index': 6,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [7.27672888504845, 4.23557550026145, 12.3922844230364],
         'tag': 0,
         'index': 7,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Ta',
         'position': [0.0, 0.0, 14.7845694537418],
         'tag': 0,
         'index': 8,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [2.93778109210234, 0.0, 14.7845694537418],
         'tag': 0,
         'index': 9,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Ta',
         'position': [1.44631589764871, 2.53102532088718, 14.7845694537418],
         'tag': 0,
         'index': 10,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'Au',
         'position': [4.38409698975104, 2.53102532088718, 14.7845694537418],
         'tag': 0,
         'index': 11,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0},
        {'symbol': 'S',
         'position': [4.99301884665953, 1.49768119401374, 16.7347550157333],
         'tag': 1,
         'index': 12,
         'charge': 0.0,
         'momentum': [0.0, 0.0, 0.0],
         'magmom': 0.0}],
       'cell': {'array': [[5.87556258420465, 0.0, 0.0],
         [2.89263179529741, 5.06205064177436, 0.0],
         [0.0, 0.0, 24.7845692460728]],
        'pbc': [True, True, True],
        '__ase_objtype__': 'cell'},
       'pbc': [True, True, True],
       'info': {},
       'constraints': [{'name': 'FixAtoms',
         'kwargs': {'indices': [0, 1, 2, 3, 4, 5, 6, 7]}}],
       'natoms': 13,
       'mass': 2299.5466939999997,
       'spacegroup': 'P1 (1)',
       'chemical_symbols': ['Au', 'S', 'Ta'],
       'symbol_counts': {'Ta': 6, 'Au': 6, 'S': 1},
       'volume': 737.1524571000787},
      'calc': {'calculator': {'module': 'ase.calculators.singlepoint',
        'class': 'SinglePointCalculator'}},
      'results': {'energy': -21550.0244970718,
       'forces': [[-0.343776547268582, 0.198238995713678, 0.91678828152204],
        [0.256798441315125, -0.139238863030386, -0.327803311768777],
        [0.254964730540496, -0.147942561451971, -0.319660885167301],
        [-0.466383715279558, 0.233086956630799, 0.972775081300305],
        [0.0579110133841712, 0.195392013183186, -0.0490314515063543],
        [-0.0205264020881259, -0.0250787973454987, 0.246483175432047],
        [0.317658510025292, -0.306177505968886, -0.286142955915514],
        [0.0754551929155217, -0.0183268233336459, 0.113915810884565],
        [-0.262860588348805, 0.418947119064149, 0.783792056967257],
        [-0.53517769412482, 0.0478441160664724, -0.126114409833286],
        [0.0452349605783995, -0.212374663784112, -0.696931193317955],
        [-1.277769879815, 1.39071587058713, -1.22322829712369],
        [1.89847172105557, -1.63508637055155, -0.00484215858364972]],
        'fmax': 1.89847172105557},
        'reconstructed': 0,
        'fw_id': 50514,
        'site': 'bridge',
        'site_type': 'A_B|A',
        'coordination': 'Au-Ta',
        'neighborcoord': ['Au:Au-Au-Au-Ta-Ta-Ta-Ta-Ta-Ta',
        'Ta:Au-Au-Au-Au-Au-Au-Ta-Ta-Ta'],
        'nextnearest': 'Au-Au-Au-Au-Au-Au-Au-Au-Ta-Ta-Ta-Ta-Ta-Ta-Ta-Ta',
        'nextnearestcoordination': 'Au-Au-Ta-Ta'},
        'coordination': 'Au-Ta-Ta',
        'neighborcoord': ['Ta:Au-Au-Au-Au-Au-Au-Ta-Ta-Ta',
        'Au:Au-Au-Au-Ta-Ta-Ta-Ta-Ta-Ta',
        'Ta:Au-Au-Au-Au-Au-Au-Ta-Ta-Ta'],
        'nextnearest': 'Au-Au-Au-Au-Au-Au-Au-Au-Au-Au-Au-Ta-Ta-Ta-Ta-Ta-Ta-Ta-Ta-Ta',
        'nextnearestcoordination': 'Au-Au-Ta-Ta'}
    
    from gplearn.site_rebuilt import SiteRebuild, SiteClassification
    from gplearn.read_datasets import make_atoms_from_doc
    iatoms = make_atoms_from_doc(doc['initial_configuration'])
    sc = SiteRebuild(iatoms,  'L10',
                 metalA =doc['metalA'], 
                 metalB =doc['metalB'],  
                 adsorbate='O')
    for site in ["A_B|B", "B_B|A", "A_A|B", "A_B_B|HCP", "A", "B", "A_A_B|FCC", \
                    "A_B_B|FCC", "A_A_B|HCP", "A_B|A"]:
        rebuild = sc.rebuild_initial(site)
        scla = SiteClassification(rebuild, A = rebuild)
        got_site = scla.get_site()[1]
        got_site = got_site.replace(doc['metalA'], "A")
        got_site = got_site.replace(doc['metalB'], "B")
        if "|" not in site:
            assert site == got_site
        else:
            _A , _B = got_site.split("|")
            __A , __B = site.split("|")

            _A = sorted(_A.split("_"))
            _B = sorted(_B.split("_"))
            __A = sorted(__A.split("_"))
            __B = sorted(__B.split("_"))
            print(_A, __A)
            assert _A == __A
            assert _B == __B
    assert _B == __B
    

def test_A1():
    doc =  {'slab_name': 'Au',
 'reconstructed': 0,
 'site': 'hollow',
 'metalA': 'Au',
 'metalB': 'nan',
 'site_type': 'A_A_A|FCC',
 'fw_id': 31066,
 'energy': -0.17403244629525716,
 'adsorbate': 'S',
 'initial_configuration': {'atoms': {'atoms': [{'symbol': 'Au',
     'position': [1.4877364947175, 0.858945039215122, 9.99999929233104],
     'tag': 0,
     'index': 0,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [4.4632094841525, 0.858945039215122, 9.99999929233104],
     'tag': 0,
     'index': 1,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [2.975472989435, 3.43578025686048, 9.99999929233104],
     'tag': 0,
     'index': 2,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [5.95094597887, 3.43578025686048, 9.99999929233104],
     'tag': 0,
     'index': 3,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [5.95094597887, 1.71789017843024, 12.4294628204054],
     'tag': 0,
     'index': 4,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [2.975472989435, 1.71789017843024, 12.4294628204054],
     'tag': 0,
     'index': 5,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [7.4386824735875, 4.2947252960756, 12.4294628204054],
     'tag': 0,
     'index': 6,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [4.4632094841525, 4.2947252960756, 12.4294628204054],
     'tag': 0,
     'index': 7,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [0.0, 0.0, 14.8589263484798],
     'tag': 0,
     'index': 8,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [2.975472989435, 0.0, 14.8589263484798],
     'tag': 0,
     'index': 9,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [1.4877364947175, 2.57683521764536, 14.8589263484798],
     'tag': 0,
     'index': 10,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'Au',
     'position': [4.4632094841525, 2.57683521764536, 14.8589263484798],
     'tag': 0,
     'index': 11,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0},
    {'symbol': 'S',
     'position': [1.4877364947175, 0.858945039215122, 16.4254262376235],
     'tag': 1,
     'index': 12,
     'charge': 0.0,
     'momentum': [0.0, 0.0, 0.0],
     'magmom': 0.0}],
   'cell': {'array': [[5.95094557887002, 0.0, 0.0],
     [2.97547278943501, 5.1536706352907, 0.0],
     [0.0, 0.0, 24.8589252408109]],
    'pbc': [True, True, True],
    '__ase_objtype__': 'cell'},
   'pbc': [True, True, True],
   'info': {},
   'constraints': [{'name': 'FixAtoms',
     'kwargs': {'indices': [0, 1, 2, 3, 4, 5, 6, 7]}}],
   'natoms': 13,
   'mass': 2395.6588279999996,
   'spacegroup': 'P1 (1)',
   'chemical_symbols': ['Au', 'S'],
   'symbol_counts': {'Au': 12, 'S': 1},
   'volume': 762.4036851443883},
  'calc': {'calculator': {'module': 'ase.calculators.singlepoint',
    'class': 'SinglePointCalculator'}},
  'results': {'energy': -19386.7146468816,
   'forces': [[-5.4507387019261e-05, -9.8987471709507e-05, 0.0017707187471776],
    [-7.71330948385768e-05, 0.00391681855590293, -0.150632192468992],
    [0.00339205640068448, -0.00176943319559695, -0.149922310886161],
    [-0.0035972304329551, -0.00200520335548687, -0.151229973953991],
    [-0.0312844119355784, 0.0180910531737559, -0.0527544088838962],
    [0.0314124528730104, 0.0182856856830653, -0.0526937308492899],
    [0.0002077451354319, 0.00012032762794818, 0.623672236722803],
    [-5.65642695482897e-06, -0.0364006501162212, -0.0525170960621096],
    [-1.30179889573975, -0.75172782944286, -0.386232145329632],
    [1.30189325522576, -0.751417240180977, -0.38614524204278],
    [-7.79044257869626e-05, 1.50293526760588, -0.386435262479373],
    [1.59408395999725e-05, 2.98247966709164e-05, -0.404608076733658],
    [-2.54539212967304e-05, 4.010920931606e-05, 1.5477274842199]],
   'fmax': 1.5477274842199},
  'reconstructed': 0,
  'fw_id': 31066,
  'site': 'hollow',
  'site_type': 'A_A_A|FCC'
                          }}
    from gplearn.site_rebuilt import SiteRebuild, SiteClassification
    from gplearn.read_datasets import make_atoms_from_doc
    iatoms = make_atoms_from_doc(doc['initial_configuration'])
    sc = SiteRebuild(iatoms,  'A1',
                 metalA =doc['metalA'], 
                 metalB =None,  
                 adsorbate='O')
    for site in ["A", "A_A_A|HCP", "A_A|A", "A_A_A|FCC"]:
        rebuild = sc.rebuild_initial(site)
        scla = SiteClassification(rebuild, A = rebuild)
        got_site = scla.get_site()[1]
        got_site = got_site.replace(doc['metalA'], "A")
        assert site == got_site
    assert site == got_site
    