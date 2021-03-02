Basic Example
=============

===
API
===

A simple example of how to use Pypka as an API is provided below. In
this snippet we are estimating the pKa values for all titrable sites
in the `4lzt.pdb structure <https://files.rcsb.org/download/4LZT.pdb>`_ downloaded from the Protein Data Bank.

.. code-block:: python
   
   from pypka import Titration
   
   params = {
    'structure'     : '4lzt.pdb',    
    'ncpus'         : -1,
    'epsin'         : 15,
    'ionicstr'      : 0.1,
    'pbc_dimensions': 0
    #Set the ffinput when using PDB files from simulations    
    #'ffinput': 'GROMOS' # options: GROMOS, AMBER, CHARMM
   }
   
   tit = Titration(params)
      
   pH = 7.0
   for site in tit:
       state = site.getProtState(pH)[0]    
       print(site.res_name, site.res_number, site.pK, state)         
   
   
You may also try it out on a `online notebook.
<https://colab.research.google.com/github/mms-fcul/PypKa/blob/master/pypka/example/notebook/pypka.ipynb>`_ 


===
CLI
===

The same calculation can be done via CLI by resorting to a input
parameter file.

.. code-block:: text
   :caption: parameter.dat
      
   structure       = 4lzt.pdb
   epsin           = 15
   ionicstr        = 0.1
   pbc_dimensions  = 0
   ncpus           = -1
   sites           = all   
   #Set the ffinput when using PDB files from simulations
   #Options: GROMOS, AMBER, CHARMM

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

   number of dimensions with periodic boundaries. 0 for solvated proteins and 2 for lipidic systems

.. object:: ncpus
	    
   :type: int

   number of CPUs to use in the calculations (-1 to use all available)
