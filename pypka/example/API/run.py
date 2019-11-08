import sys
sys.path.insert(1, '../../')

from pypka import Titration

parameters = {'structure'     : 'lyso.pdb',
              'convergence'   : 0.01,
              'pH'            : "0,14",
              'scaleM'        : 2,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : sys.argv[1],
              'output'        : 'pKas.out',
              'titration_output': 'titration.out',
              }
              #'gsize'         : 91,
              #'bndcon'        : 4,
              #'nanoshaper'    : 0,
              #'precision'     : 'single'}
pKa = Titration(parameters)

#for pH in pKa.getTitrationCurve('total').keys():
#    print(pH, pKa.getTitrationCurve('total')[pH])
exit()
for site in pKa:
    print(site, pKa[site], pKa.getProtState(site, 7))
    print(pKa.getMostProbStates(site, 7))
    print(pKa.getStatesProb(site, 7))
    print(pKa.getFinalState(site, 7))
    print(pKa.getTitrationCurve(site))
    exit()
