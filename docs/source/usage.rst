Usage
===========

.. toctree::
   :maxdepth: 2
   
   usage/methods
   usage/parameters


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
