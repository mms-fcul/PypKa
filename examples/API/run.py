from pypka import Titration

params = {
    "structure": "4lzt.pdb",  # Input structure file name
    "ncpus": -1,  # Number of processes (-1 for maximum allowed)
    "epsin": 15,  # Dielectric constant of the protein
    "ionicstr": 0.1,  # Ionic strength of the medium (M)
    "pbc_dimensions": 0,  # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
}

tit = Titration(params)

pH = 7.0
for site in tit:
    state = site.getProtState(pH)[0]
    print(site.res_name, site.res_number, site.pK, state)
