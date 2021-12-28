Using a MD structure
========================

In this example we are calculating the p\ :emphasis:`K`\ \ :sub:`a`\  values for all titrable residues using a molecular dynamics generated structure.
The argument :code:`ffinput` needs to be set to the appropriate nomenclature scheme.


.. code-block:: text
   :caption: parameters.dat
      
   structure       = amber.pdb     # Input structure file name
   ncpus           = -1            # Number of processes (-1 for maximum allowed)
   epsin           = 15            # Dielectric constant of the protein
   ionicstr        = 0.1           # Ionic strength of the medium (M)
   pbc_dimensions  = 0             # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
   sites           = all           # Titrate all available sites
   output          = pKas.out      # pKa values output file
   ffinput         = AMBER         # Input structure nomenclature scheme ("AMBER", "CHARMM", "GROMOS")

To run the simulation execute the following command:

.. code-block:: bash

   python3 -m pypka parameters.dat

