p\ :emphasis:`K`\ \ :sub:`a`\  Calculation
==========================================

In this example we are estimating the p\ :emphasis:`K`\ \ :sub:`a`\  values for all titrable sites of a hen egg-white lysozyme.

The structure may be downloaded from the Protein Data Bank::

   wget https://files.rcsb.org/download/4lzt.pdb


PypKa features a API and a CLI. Both interfaces provide the same amount of features.

API
---

.. code-block:: python
   
   from pypka import Titration
   
   params = {
    'structure'     : '4lzt.pdb', # Input structure file name
    'ncpus'         : -1,         # Number of processes (-1 for maximum allowed)
    'epsin'         : 15,         # Dielectric constant of the protein
    'ionicstr'      : 0.1,        # Ionic strength of the medium (M)
    'pbc_dimensions': 0           # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
   }
   
   tit = Titration(params)
      
   pH = 7.0
   for site in tit:
       state = site.getProtState(pH)[0]    
       print(site.res_name, site.res_number, site.pK, state)   
   
   
You may also try it out on a `online notebook.
<https://colab.research.google.com/github/mms-fcul/PypKa/blob/master/pypka/example/notebook/pypka.ipynb>`_ 



CLI
---

The same calculation can be done via CLI by resorting to a input
parameter file.

.. code-block:: text
   :caption: parameters.dat
      
   structure       = 4lzt.pdb      # Input structure file name
   ncpus           = -1            # Number of processes (-1 for maximum allowed)
   epsin           = 15            # Dielectric constant of the protein
   ionicstr        = 0.1           # Ionic strength of the medium (M)
   pbc_dimensions  = 0             # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
   sites           = all           # Titrate all available sites
   output          = pKas.out      # pKa values output file

To run the simulation execute the following command:

.. code-block:: bash

   python3 -m pypka parameters.dat

