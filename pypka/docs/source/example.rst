Basic Example
=============

===
API
===

A simple example of how to use Pypka as an API is provided below. In
this snippet we are estimating the pKa values for all titrable sites
in the 4lzt.pdb structure.

.. code-block:: python
   
   from pypka.pypka import Titration

   parameters = {'structure'     : '4lzt.pdb',
                 'epsin'         : 10,
                 'ionicstr'      : 0.1,
                 'pbc_dimensions': 0,
                 'ncpus'         : 4}
   
   pKa = Titration(parameters)
   
   for site in pKa:
       print(site, pKa[site], pKa.getProtState(site, 7))
   
   print(pKa.getParameters())


===
CLI
===

The same calculation can be done via CLI by resorting to a input
parameter file.

.. code-block:: text
   :caption: parameter.dat
      
   structure       = 4lzt.pdb
   epsin           = 10
   ionicstr        = 0.1
   pbc_dimensions  = 0
   ncpus           = 4
   sites_A         = all

To execute pypka simply type one of the two:

.. code-block:: bash

   pypka parameters.dat
   python3 -m pypka parameters.dat

====================
Mandatory Parameters
====================

.. object:: structure
	    
   :type: str

   the PDB filename

.. object:: epsin
	    
   :type: float

   internal dielectric to be used in the PB calculations.

.. object:: ionicstr
	    
   :type: float

   ionic strength of the medium

.. object:: pbc_dimensions
	    
   :type: int

   number of dimensions with periodic boundaries. 0 for solvated
   proteins and 2 for lipidic systems

.. object:: ncpus
	    
   :type: int

   number of CPUs to use in the calculations
