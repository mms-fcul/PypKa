import sys
sys.path.insert(1, '../../')

from pypka import Titration

parameters = {'structure'     : 'lyso.pdb',
              'gsize'         : 41,
              'convergence'   : 0.2,
              'pH'            : "6,8",
              'scaleM'        : 2,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : 2,
              'output'        : 'pKas.out',
              'titration_output': 'titration.out',
              }
              #'bndcon'        : 4,
              #'nanoshaper'    : 0,
              #'precision'     : 'single'}
pKa = Titration(parameters)

for site in pKa:
    print(site, pKa[site], pKa.getProtState(site, 7))
    print(pKa.getMostProbStates(site, 7))
    print(pKa.getStatesProb(site, 7))
    print(pKa.getFinalState(site, 7))
    print(pKa.getTitrationCurve(site))
    exit()