import sys
sys.path.insert(1, '../../')

from pypka import Titration

parameters = {'structure'     : 'lyso.pdb',      
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : 2,
              'output'        : 'pKas.out'}
pKa = Titration(parameters)
