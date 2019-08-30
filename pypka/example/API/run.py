import sys
sys.path.insert(1, '../../')

from pypka import Titration

parameters = {'structure'     : 'lyso.pdb',
              'gsize'         : 81,
              'scaleM'        : 1,
              'nanoshaper'    : 0,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'bndcon'        : 4,
              'ncpus'         : 2,
              'output'        : 'pKas.out',
              'debug'         : True,
              'precision'     : 'double'}
pKa = Titration(parameters)

