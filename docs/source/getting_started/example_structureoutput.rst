MD Structure Preparation
========================

In this example we are preparing a structure for molecular dynamics by assigning the most likely protonation state to all titrable sites.

In order to output this structure we need to add the argument :code:`structure_output`. Since we aren't interested in calculating p\ :emphasis:`K`\ \ :sub:`a`\  values, we can speed up the simulation by setting :code:`ph` to the single pH value of the desired output structure.


.. code-block:: text
   :caption: parameters.dat
      
   structure       = 4lzt.pdb      # Input structure file name
   ncpus           = -1            # Number of processes (-1 for maximum allowed)
   epsin           = 15            # Dielectric constant of the protein
   ionicstr        = 0.1           # Ionic strength of the medium (M)
   pbc_dimensions  = 0             # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
   sites           = all           # Titrate all available sites
   output          = pKas.out      # pKa values output file
   structure_output = out.pdb, 7, amber  # Output structure with the most likely protonation state at a given pH value
   # filename, pH value, force field ("gromos_cph", "amber")
   pH = 7                          # pH range

To run the simulation execute the following command:

.. code-block:: bash

   python3 -m pypka parameters.dat

