# Adding extra residues

For demonstration purposes, we will add an FMOC residue to the CHARMM36m FF for which the parameters have been kindly supplied by Alexander van Teijlingen.

DataBaseT.crg and DataBaseT.siz are the main files where the charges and radii of atoms have to be added.

The new DataBaseT files will be generated by extend_db.py (tmp.crg and tmp.siz) provided one follows the steps describe below.

If you believe your residue could benefit other users please create a pull request on GitHub or send the new parameters to the developers directly.


### Add residue block to ffG54a7pHt.rtp file

The lines should adhere to GROMACS .rtp file logic. Please use unique atom types to describe your residues in order to avoid clashes with previously declared atoms.

```
# EXTRA RESIDUES

[ CHO ] ; - adapted from ATB entry M81X
 [ atoms ]
     H3'  H       0.450 0
      O3  OAlc   -0.710 0
      C3  C       0.230 0
      H3  HC      0.040 0
      C1  CH2     0.080 0
      C2  CH2    -0.080 0
      C4  CH2    -0.010 0
```


### Add atom types to ffG54a7pHtnb.itp

Any new atom types must also be added to the nonbonded parameters file under the section "nonbond_params". Only the interaction with OW is necessary, in order to calculate radii.

```
  OW     OAlc 1  0.002155371     1.7853e-06      ; ATB entry M81X
```


### Add placeholder info in DataBaseT_old.crg and DataBaseT_old.siz

```
H3'   CHO       X.XXX
O3    CHO       X.XXX
C3    CHO       X.XXX
H3    CHO       X.XXX
C1    CHO       X.XXX
C2    CHO       X.XXX
C4    CHO       X.XXX
```


### Run extend_db

```
python3 extend_db.py
```

New tmp.crg and tmp.siz have been generated. Please inspect them to check the new residue has been added at the bottom correctly. After the visual inspection, rename them DataBaseT.crg and DataBaseT.siz, respectively.

```
mv tmp.crg DataBaseT.crg
mv tmp.siz DataBaseT.siz
```


### Add residue to pdb2pqr

The residue also needs to be added to PDB2PQR's database of residues, so that is not discarded on the preprocessing step.

/pdb2pqr/dat/CHARMM.DAT
```

# FMOC - Alexander van Teijlingen 07-12-2020
FMO     C1    0.01  0.01    C1
FMO     H1    0.01  0.01    H1
FMO     C2    0.01  0.01    C2
FMO     H2    0.01  0.01    H2
```