import readFiles
import rundelphi

from ctypes import *

import os


# if using parallel version don't forget to set system-wide variables
# export OMP_NUM_THREADS=8
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pedror/delphit/dependencies/NanoShaper0.7/build_lib:/home/pedror/delphit/dependencies/cppDelphi77/lib/

# Precision
precision = 'single'

# PASS ARGUMENTS

# use perfil or igrid
# set the other to 0
py_perfil  = 0
py_igrid   = 121
py_scale   = 1.0223857030
py_scale1  = py_scale    
py_repsin  = 2.0
py_repsout = 80.0
py_radprb  = 1.4
py_conc    = 0.1
py_ibctyp  = 4
py_res2    = 0.01
py_nlit    = 50

py_acent = [44.014700,44.014700,114.781000]

# custom.prm and surfaceConfiguration.prm have to exist
# nanoshaper can not be used with focussing
# because the grid has to be bigger than the system
# -1 to turn off Nanoshaper
# 0 for connoly surface
py_isurftype = -1

# cpp new solver
py_parallel = False

# single characters only for energy and site arguments

py_energy = ['s', 'c']

py_site = ['a', 'q', 'p']

py_in_crg = "P.crg"
py_in_crg_len = len(py_in_crg)

py_in_siz = "DataBaseT.siz"
py_in_siz_len = len(py_in_siz)

py_in_pdb = "P.pdb"
py_in_pdb_len = len(py_in_pdb)

natom = 17250

if precision == 'double':
    float_type = c_double
else:
    float_type = c_float
    
py_atpos   = float_type * 3 * natom
p_atpos    = py_atpos()
py_i_atpos = addressof(p_atpos)

py_rad3   = float_type * natom
p_rad3    = py_rad3()
py_i_rad3 = addressof(p_rad3)

py_chrgv4   = float_type * natom
p_chrgv4    = py_chrgv4()
py_i_chrgv4 = addressof(p_chrgv4)

py_atinf   = (c_char * 15  * natom)()
py_i_atinf = addressof(py_atinf)

nmedia = 1
nobject = 1
len_medeps = nmedia + nobject
py_medeps  = float_type * len_medeps
p_medeps    = py_medeps()
py_i_medeps = addressof(p_medeps)

len_iatmed  = natom + 1
py_iatmed   = c_int * len_iatmed
p_iatmed    = py_iatmed()
py_i_iatmmed = addressof(p_iatmed)

py_dataobject   = (c_char * 96 * nobject * 2)()
py_i_dataobject = addressof(py_dataobject)

py_rmaxdim = 999.999

py_rmaxdim = readFiles.delphi(py_igrid, py_scale, py_repsin, py_repsout,
                              py_acent, py_in_pdb, py_in_crg, py_in_siz,
                              natom, nobject, py_i_atpos, py_i_rad3,
                              py_i_chrgv4, py_i_atinf, py_i_medeps,
                              py_i_iatmmed, py_i_dataobject, py_rmaxdim)

#print natom
#for i in range(natom):
#    print p_atpos[i][0], p_atpos[i][1], p_atpos[i][2] #, p_sitpot[i]

if py_perfil != 0:
    py_igrid = int(py_scale * 100 / py_perfil * py_rmaxdim)
    if py_igrid % 2 == 0:
        py_igrid += 1

#py_nonit   = 5
#py_relfac  = 0.200000
#py_relpar  = 0.750000
#py_pbx     = True
#py_pby     = True

py_nonit   = 0
py_relfac  = 0.0
py_relpar  = 0.0
py_pbx     = False
py_pby     = False

py_in_frc  = 'self'

esolvation = 999.999

py_sitpot   = float_type * natom
p_sitpot    = py_sitpot()
py_i_sitpot = addressof(p_sitpot)


py_out_phi  = True

len_phimap   = py_igrid * py_igrid * py_igrid
py_phimap4   = c_float * len_phimap
p_phimap4    = py_phimap4()
py_i_phimap4 = addressof(p_phimap4)

esolvation = rundelphi.delphi(py_igrid, py_scale, py_repsin, py_repsout,
                              py_radprb, py_conc, py_ibctyp, py_res2,
                              py_nlit, py_acent, py_energy, py_site,
                              py_nonit, py_relfac, py_relpar, py_pbx,
                              py_pby, py_in_frc, natom, nmedia, nobject,
                              py_i_atpos, py_i_rad3, py_i_chrgv4,
                              py_i_atinf, py_i_medeps, py_i_iatmmed,
                              py_i_dataobject, py_i_phimap4, py_scale1,
                              py_out_phi, py_i_sitpot, esolvation,
                              py_isurftype, py_parallel)

#c = 0
#ind=-1
#for i in range(py_igrid):
#    for ii in range(py_igrid):
#        for iii in range(py_igrid):
#            ind += 1
#            if p_phimap4[ind] != 0.0:
#                c+=1
#                if c < 20:
#                    print p_phimap4[ind]
print(esolvation)


py_scale   = 4.0895400000
py_ibctyp  = 3
py_nlit    = 500

py_out_phi = False

py_nonit   = 0
py_relfac  = 0.0
py_relpar  = 0.0
py_pby     = False
py_pbx     = False


esolvation = rundelphi.delphi(py_igrid, py_scale, py_repsin, py_repsout,
                              py_radprb, py_conc, py_ibctyp, py_res2,
                              py_nlit, py_acent, py_energy, py_site,
                              py_nonit, py_relfac, py_relpar, py_pbx,
                              py_pby, py_in_frc, natom, nmedia, nobject,
                              py_i_atpos, py_i_rad3, py_i_chrgv4,
                              py_i_atinf, py_i_medeps, py_i_iatmmed,
                              py_i_dataobject, py_i_phimap4, py_scale1,
                              py_out_phi, py_i_sitpot, esolvation,
                              py_isurftype, py_parallel)

#for i in range(natom):
#    if p_sitpot[i] != 0.0:
#        print p_atpos[i][0], p_atpos[i][1], p_atpos[i][2], p_sitpot[i]

#for i in range(natom):
#    print py_atinf[i].value

#print p_medeps[:]
#print p_iatmed[:]

#for i in range(nobject):
#    print py_dataobject[i][0].value

print(esolvation)

print('exiting')
