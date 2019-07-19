import readFiles
from ctypes import *

# Precision
precision='single'

# PASS ARGUMENTS
py_igrid   = 81
py_scale   = 4.0
py_repsin  = 2.0
py_repsout = 80.0
py_radprb  = 1.4
py_conc    = 0.1
py_ibctyp  = 4
py_res2    = 0.01
py_nlit    = 500

py_acent = [22.855, 27.190, 9.5600]

# single characters only for energy and site arguments

py_energy = ['s', 'c']

#py_site = ['a', 'q', 'p']
py_site = []

py_in_crg = "DataBaseT.crg"
py_in_siz = "DataBaseT.siz"
py_in_pdb = "1001-Tyr1N.pdb"

natom = 14

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

py_rmaxdim = readFiles.delphi(py_igrid, py_scale, py_repsin,
                              py_repsout, py_acent, py_in_pdb,
                              py_in_crg, py_in_siz, natom, nobject,
                              py_i_atpos, py_i_rad3, py_i_chrgv4,
                              py_i_atinf, py_i_medeps, py_i_iatmmed,
                              py_i_dataobject, py_rmaxdim)
print natom
for i in range(natom):
    print p_atpos[i][0], p_atpos[i][1], p_atpos[i][2]

#, p_rad3[i], p_chrgv4[i]    
for i in range(natom):
    print py_atinf[i].value

print p_medeps[:]
print p_iatmed[:]

for i in range(nobject):
    print py_dataobject[i][0].value

print py_rmaxdim
print 'exiting'
