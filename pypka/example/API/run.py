import sys

sys.path.insert(1, "../../")

import pypka

# from pypka import Titration

parameters = {
    "structure": "lyso.pdb",
    "convergence": 0.1,
    "pH": "0,14",
    "scaleM": 2,
    "epsin": 10,
    "ionicstr": 0.1,
    "pbc_dimensions": 0,
    "ncpus": sys.argv[1],
    "output": "pKas.out",
    "titration_output": "titration.out",
    "ser_thr_titration": False,
    "structure_output": ("out4.pdb", 4.5, "gromos_cph"),
}
#'gsize'         : 91,
#'bndcon'        : 4,
#'nanoshaper'    : 0,
#'precision'     : 'single'}

sites = {"A": ["1N", 1, 7, 15, 18, 35, 48, 52, 66, 87, 101, 119, "129C"]}

pKa = pypka.Titration(parameters, sites=sites)
# exit()

# Run each step separately
# pKa = Titration(parameters, run='preprocess')
# pKa.DelPhiLaunch()
# pKa.calcSiteInteractionsParallel()
# pKa.run_mc()
# print(pKa)

# from IPython import embed; embed()

# for pH in pKa.getTitrationCurve('total').keys():
#    print(pH, pKa.getTitrationCurve('total')[pH])

# for site in pKa:
#    print(site, pKa[site], pKa.getProtState(site, 7))
#    print(pKa.getMostProbStates(site, 7))
#    print(pKa.getStatesProb(site, 7))
#    print(pKa.getFinalState(site, 7))
#    print(pKa.getTitrationCurve(site))
