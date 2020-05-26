# DOWNLOAD PDB files
#from Bio.PDB import PDBList
#import os
#pdbl = PDBList()
#pdbl.retrieve_pdb_file('1LSE', file_format="pdb", pdir=".")
#pdbl.retrieve_pdb_file('4LZT', file_format="pdb", pdir=".")
#os.rename('pdb1lse.ent', '1lse.pdb')
#os.rename('pdb4lzt.ent', '4lzt.pdb')

import sys
sys.path.insert(1, '../../')
from pypka import Titration

parameters = {'structure'     : '4lzt.pdb',      
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'temp'          : 100,
              'grid_fill'     : 0.8,         # FUTURE VERSION
              'ncpus'         : 4,
              'pH'            : '0, 15',
              'maxc'          : 0.01,
              'pHstep'        : 0.5,
              'output': 'pkas.out'}

sites = {'A': ('1N', '18', '35', '48', '66', '129C')}
pKa = Titration(parameters, sites=sites, debug=True)

for site in pKa:
    print site, pKa[site], pKa.getProtState(site, 7)

print pKa.getParameters()
