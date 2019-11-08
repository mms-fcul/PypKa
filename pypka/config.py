from copy import copy
import os

tit_mole = None
sites    = {}
pid      = None
debug    = False

input_conversion = {'grid_fill': 'perfil',
                    'pbc_dimensions': 'pbc_dim',
                    'convergence': 'maxc'}
stdout = None
stderr = None

stdout_file = None
stderr_file = None

njobs   = None # to be redeclared in concurrency
pb_time = None # to be redeclared in concurrency
total_jobs = None # to be redeclared in concurrency

# DelPhi related
params = {}
default_params = {'perfil': 0.9,
                  'gsize': 81,
                  'scaleP': 1,
                  'scaleM': 4,
                  'precision': 'single',
                  'ionicstr': 0.1,
                  'bndcon': 3,
                  'maxc': 0.01,
                  'nlit': 500,
                  'nonit': 0,
                  'relfac': 0.75,
                  'relpar': 0.75,
                  'pbx': False,
                  'pby': False,
                  'ffID': 'G54A7',
                  'ffinput': 'GROMOS',
                  'pHmin': 0,
                  'pHmax': 14,
                  'pHstep': 0.25,
                  'epssol': 80.0,
                  'temp': 298.0,
                  'seed': 1234567,
                  'cutoff': 2.5,
                  'pbc_dim': 0,
                  'epsin': 20.0,
                  'slice': 0.05,
                  'ncpus': 1,
                  'clean_pdb': True,
                  'keep_ions': False,
                  'couple_min': 2.0,
                  'mcsteps': 200000,
                  'eqsteps': 1000}

# -1 NanoShaper off
# 0 connolly surface
# 1 skin
# 2 blobby
# 3 mesh
# 4 msms
nanoshaper = -1

# Paths
# TODO: include in dependencies
# in dependencies there is a pdb2pqr that has a important dat folder
fileDir = os.path.dirname(os.path.abspath(__file__))
script_dir = fileDir
pdb2pqr = "{0}/pdb2pqr/pdb2pqr.py".format(fileDir)
userff = "{0}/pdb2pqr/dat/GROMOS.DAT".format(fileDir)
usernames = "{0}/pdb2pqr/dat/GROMOS.names".format(fileDir)


# Input Files
f_in = None
f_in_extension = None
f_out = None
f_prot_out = None
f_structure_out = None
f_log = "LOG"
f_dat = None

# Force Field Files
f_crg = None
f_siz = None

lipids = {'cholesterol': 'CHO',  # to edit
          'POPC': 'POP'}
lipid_residues = ['DMX', 'POX', 'PJ2', 'PJ1', 'PJ0', 'CHL']  # allowed residue names

ions = ['CL-', 'NA+']


# Constants
kBoltz = 5.98435e-6  # e^2/(Angstrom*K)
log10 = 2.302585092994046


# Tautomer Variables
terminal_offset = 2000

TITRABLETAUTOMERS = {'LYS': 3,
                     'HIS': 2,
                     'ASP': 4,
                     'GLU': 4,
                     'SER': 3,
                     'THR': 3,
                     'CYS': 3,
                     'CTR': 4,
                     'NTR': 3,
                     'TYR': 2}


TITRABLERESIDUES = list(TITRABLETAUTOMERS.keys())
REGULARTITRATINGRES = copy(list(TITRABLETAUTOMERS.keys()))

for res in REGULARTITRATINGRES:
    ntautomers = TITRABLETAUTOMERS[res]
    for i in range(ntautomers):
        TITRABLERESIDUES.append(res[0:2] + str(i + 1))


gromos2amber = {
    'ASP': {0: {'HD11': 'HD2', 
                'OD1': 'OD2', 
                'OD2': 'OD1'},
            1: {'HD21': 'HD2'},
            2: {'HD12': 'HD2', 
                'OD1': 'OD2', 
                'OD2': 'OD1'},
            3: {'HD22': 'HD2'}
    },
    'CYS': {0: {'HG1': 'HG'},
            1: {'HG2': 'HG'},
            2: {'HG3': 'HG'}
    },
    'GLU': {0: {'HE11': 'HE2',
                'OE1': 'OE2', 
                'OE2': 'OE1'},
            1: {'HE21': 'HE2'},
            2: {'HE12': 'HE2',
                'OE1': 'OE2', 
                'OE2': 'OE1'},
            3: {'HE22': 'HE2'}
    },
    'HIS': {},
    'TYR': {0: {'HH1': 'HH'},
            1: {'HH2': 'HH'}
    },
    'LYS': {0: {'HZ3': 'HZ1'},
            1: {'HZ3': 'HZ2'}},
    'SER': {0: {'HG1': 'HG'},
            1: {'HG2': 'HG'},
            2: {'HG3': 'HG'}
    },
    'THR': {0: {'HG1': 'HG1'},
            1: {'HG2': 'HG1'},
            2: {'HG3': 'HG1'}
    },
    'NTR': {0: {'H3': 'H1'},
            1: {'H3': 'H2'}},
    'CTR': {0: {'HO11': 'HO',
                'O1': 'O',
                'O2': 'OXT',
                'CT': 'C'},
            1: {'HO21': 'HO',
                'O2': 'O',
                'O1': 'OXT',
                'CT': 'C'},
            2: {'HO12': 'HO',
                'O1': 'O',
                'O2': 'OXT',
                'CT': 'C'},
            3: {'HO22': 'HO',
                'O2': 'O',
                'O1': 'OXT',
                'CT': 'C'},
            4: {'O1': 'O',
                'O2': 'OXT',
                'CT': 'C'}
    }
}


ffconversions = {'GROMOS': {'AMBER': gromos2amber}}

AMBER_Hs = {
    'NTR': ('H1', 'H2', 'H3'),
    'CTR': ('HO'),
    'ASP': ('HD2'),
    'CYS': ('HG'),
    'GLH': ('HE2'),
    'HIP': ('HD1', 'HE2'),
    'HID': ('HD1'),
    'HIE': ('HE2'),
    'LYS': ('HZ1', 'HZ2', 'HZ3'),
    'LYN': ('HZ2', 'HZ3'),
    'SER': ('HG'),
    'THR': ('HG1'),
    'TYR': ('HH')
}

AMBER_mainchain_Hs = ['H', 'HA']
mainchain_Hs = {}

AMBER_protomers = {'ASP': {'ASH': {0: ('HD21', 'HD12', 'HD22'), 1: ('HD11', 'HD12', 'HD22'), 
                                   2: ('HD11', 'HD21', 'HD22'), 3: ('HD11', 'HD21', 'HD21')}, 
                           'ASP': {4: ('HD11', 'HD12', 'HD21', 'HD22')}},
                   'CYS': {'CYS': {0: ('HG2', 'HG3'), 1: ('HG1', 'HG3'), 2: ('HG1', 'HG2')}, 
                           'CYM': {3: ('HG1', 'HG2', 'HG3')}},
                   'GLU': {'GLH': {0: ('HE21', 'HE12', 'HE22'), 1: ('HE11', 'HE12', 'HE22'), 
                                   2: ('HE11', 'HE21', 'HE22'), 3: ('HE11', 'HE12', 'HE21')}, 
                           'GLU': {4: ('HE11', 'HE12', 'HE21', 'HE22')}},
                   'HIS': {'HID': {0: ('HE2')}, 
                           'HIE': {1: ('HD1')}, 
                           'HIP': {2: ('')}}, 
                   'TYR': {'TYR': {0: ('HH2'), 1:('HH1')}, 
                           'TYM': {2: ('HH1', 'HH2')}},
                   'LYS': {'LYN': {0: ('HZ1'), 1: ('HZ2'), 2: ('HZ3')}, 
                           'LYS': {3: ('')}}, 
                   'SER': {'SER': {0: ('HG2', 'HG3'), 1: ('HG1', 'HG3'), 2: ('HG1', 'HG2')}, 
                           'SEM': {3: ('HG1', 'HG2', 'HG3')}},
                   'THR': {'THR': {0: ('HG2', 'HG3'), 1: ('HG1', 'HG3'), 2: ('HG1', 'HG2')}, 
                           'THM': {3: ('HG1', 'HG2', 'HG3')}},
                   'NTR': {'NTR': {0: ('H1'), 1: ('H2'), 2: ('H3')},
                           'NTN': {3: ('')}},
                   'CTR': {'CTH': {0: ('HO21', 'HO12', 'HO22'), 1: ('HO11', 'HO12', 'HO22'), 
                                   2: ('HO11', 'HO21', 'HO22'), 3: ('HO11', 'HO12', 'HO21')},
                           'CTR': {4: ('HO11', 'HO12', 'HO21', 'HO22')}},
}
