[![GitHub version](https://badge.fury.io/gh/mms-fcul%2Fpypka.svg)](https://badge.fury.io/gh/mms-fcul%2Fpypka) [![CircleCI](https://circleci.com/gh/mms-fcul/PypKa.svg?style=svg)](https://circleci.com/gh/mms-fcul/PypKa) [![Codacy Badge](https://api.codacy.com/project/badge/Coverage/77db3bc226c94625acd3cea0e14c23ad)](https://www.codacy.com/app/pedrishi/PypKa?utm_source=github.com&utm_medium=referral&utm_content=mms-fcul/PypKa&utm_campaign=Badge_Coverage) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/77db3bc226c94625acd3cea0e14c23ad)](https://www.codacy.com/app/pedrishi/PypKa?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mms-fcul/PypKa&amp;utm_campaign=Badge_Grade) [![Documentation Status](https://readthedocs.org/projects/pypka/badge/?version=latest)](https://pypka.readthedocs.io/en/latest/?badge=latest)


# PypKa

A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism


# Dependencies 

  - libgfortran4
  - gawk 
  - pytest
  - numpy


# License

  pypka is distributed under a LGPL-3.0, however delphi4py depends on
  DelPhi which is proprietary. To use DelPhi the user is required to
  download the DelPhi license
  [here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

# Documentation

  Documentation can be found [here](https://pypka.readthedocs.io/en/latest/). (Under development)

# Installation

  pip3 install pypka

  apt install gawk gcc gfortran libgfortran4

# Contacts

  Please submit a github issue to report bugs and to request new features.
  Alternatively you may find the developer [here](mailto:pdreis@fc.ul.pt). Please visit ou [website](http://mms.rd.ciencias.ulisboa.pt/) for more information.


# Change Log
## v0.7
  - bug fixes

## v0.6
  - python3 porting
  - dna support
  - nanoshaper integration

## v0.5
  - MC parallelization
  - delphi4py distributed within

## v0.4
  - API development
  - code level documentation improvement
  - documentation wiki developement with tutorial
  - test suite implementation

## v0.3
  - integration with delphi4py
  - integration with PDB2PQR
  - additional testing

## v0.2
  - DelPhi is now divided into two python modules imported by delphiT.py
  - input arguments are defined in delphiT.py
  
  - DelPhi routines for reading input files are bundled in delphi1_module
  - delphi1_module returns DelPhi's internal data structure

  - delphi2_module does not do any I/O
  - delphi2_module returns solvation energies and potential map

  - delphiT.py calculates intrinsic pka and interactions energies
  - delphiT.py supports multiprocessing and Nanoshaper usage
  - delphiT.py supports membrane systems
  - delphiT.py supports double/single precision

## v0.1
  - migration from DelPhi 5.0 to DelPhi 5.1_Patched
  - delphiT calls delphiT.py
  - delphiT.py calls DelPhi 5.1_Patched
  - added benchmark process and small test set

## original version
  - delphiT calls directly DelPhi 5.0
