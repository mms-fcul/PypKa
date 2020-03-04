import sys
sys.path.insert(1, '../../')
from pypka import Titration

pH = 7
parameters = {'structure'     : 'lyso_cleaned.pdb',
              'pH'            : "0,12",
              'scaleM'        : 2,
              'epsin'         : 10,
              'ionicstr'      : 0.1,
              'pbc_dimensions': 0,
              'ncpus'         : 2,
              'output'        : 'pKas.out',
              'titration_output': 'titration.out',
              'structure_output': ('gromos.pdb', pH, 'amber'),
              'gsize'         : 41,
              'convergence'   : 0.2,
              'pHstep'        : 0.5,
              'clean_pdb': False
}

#pKa = Titration(parameters, sites={'A': ['1N', '1', '7', '129C']})


tit = Titration(parameters)

#for site in tit:
#    pK = site.pK
#    taut = site.getMostProbTaut(pH)
#    state = site.getProtState(pH)
#
#    print(site.molecule.chain, site.res_name, site.res_number, pK, taut, state)
#
#    print(site.getFinalState(pH))
#    print(site.getTitrationCurve())

print(tit.getTitrationCurve())
print(tit.getParameters())
