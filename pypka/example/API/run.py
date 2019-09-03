import sys
sys.path.insert(1, '../../')

from pypka import Titration

parameters = {'structure'     : 'lyso.pdb',
              'gsize'         : 81,
              'scaleM'        : 1,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : 2,
              'output'        : 'pKas.out'}
              #'bndcon'        : 4,
              #'nanoshaper'    : 0,
              #'precision'     : 'single'}
pKa = Titration(parameters)

