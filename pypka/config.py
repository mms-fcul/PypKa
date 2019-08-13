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
params = {'perfil': 0.9,
          'gsize': 81,
          'scaleP': 1,
          'scaleM': 4,
          'precision': 'single',
          'ionicstr': 0.1,
          'bndcon': 4,
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
          'couple_min': 2.0,
          'mcsteps': 200000,
          'eqsteps': 1000}

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
f_log = "LOG"
f_dat = None

# Force Field Files
f_crg = None
f_siz = None

lipids = {'cholesterol': 'CHO',  # to edit
          'POPC': 'POP'}
lipid_residues = ['POX', 'PJ2', 'CHL']  # allowed residue names

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
