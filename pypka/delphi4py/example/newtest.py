import sys

sys.path.insert(0, "../")

from delphi4py import DelPhi4py as Delphi

# import os

f_crg = "DataBaseT.crg"
f_siz = "DataBaseT.siz"
fpdb = "P.pdb"

delphimol = Delphi(
    f_crg,
    f_siz,
    fpdb,
    121,
    4.0895400000,
    "single",
    conc=0.1,
    ibctyp=4,
    res2=0.01,
    nlit=500,
    pbx=True,
    pby=True,
    outputfile="LOG_readFiles",
)

natoms = delphimol.natoms
p_atpos = delphimol.p_atpos

delphimol.runDelPhi(
    scale_prefocus=1,
    scale=4,
    nlit_prefocus=500,
    nonit=50,
    nlit=500,
    acent=[56.21, 111.41, 5.36],
    nonit_focus=0,
    relfac_focus=0.0,
    relpar_focus=0.0,
    relpar=0.2,
    relfac=0.2,
    pbx_focus=False,
    pby_focus=False,
    debug=True,
    outputfile="LOG_runDelPhi",
)

print((delphimol.getSolvation()))  # float
# print delphimol.getSitePotential() # array

print("exiting")
