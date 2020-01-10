import sys
sys.path.insert(1, '../../')
from pypka import Titration

pH = 1
parameters = {'structure'     : 'lyso_cleaned.pdb',
              'pH'            : "0,12",
              'scaleM'        : 2,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : 16,
              'output'        : 'pKas.out',
              'titration_output': 'titration.out',
              'structure_output': ('gromos.pdb', pH),
              'gsize'         : 41,
              'convergence'   : 0.2,
              'pHstep'        : 0.5,
              'clean_pdb': False
}

#pKa = Titration(parameters, sites={'A': ['1N', '1', '7', '129C']})


pKa = Titration(parameters)

for site in pKa:
    pK = pKa[site]
    if pK != '-':
        taut = pKa.getMostProbState(site, pH)
        print(site, pKa[site], pKa.getProtState(site, pH),
              taut, pKa.getStatesProb(site, pH),
              pKa.getStateProb(site, taut, pH))
    # print(pKa.getFinalState(site, pH))
    # print(pKa.getTitrationCurve(site))
