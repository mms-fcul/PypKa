Complete Parameters File Template
=================================

.. code-block:: text
   :caption: parameters.dat
      
   # Mandatory
   structure = 4lzt.pdb    # Input structure file name
   epsin = 15              # Dielectric constant of the protein
   ionicstr = 0.1          # Ionic strength of the medium (M)
   pbc_dimensions = 0      # PB periodic boundary conditions (0 for solvated proteins and 2 for lipidic systems)
   ncpus = -1              # Number of processes (-1 for maximum allowed)
   sites = all             # Titrate all available sites
   
   # Output
   output = pKas.txt                       # pKa values output file
   titration_output = titration.txt        # Titration curves output file
   titration_curve = titration_curve.txt   # Full titration curve (whole system charge) output file
   isoelectric_point = yes                 # Calculate isoelectric point
   save_pdb = pb_input.pdb                 # save PB input file
   structure_output = out.pdb, 7, amber    # Output structure with the most likely protonation state at a given pH value
   # filename, pH value, force field ("gromos_cph", "amber")
   
   
   # Force Field
   # f_crg = /home/user/database.crg     # DelPhi charges database file path
   # f_siz = /home/user/database.siz     # DelPhi radii database file path
   # path to search for .st files -> $ffs_dir/$ffID/$/sts
   # ffs_dir = /home/user/database # default value is an internal folder
   # ffID = G54A7 # built-in options: G54A7, CHARMM36m
   # sts = sts
   
   
   # Preprocessing
   ffinput = GROMOS          # Input structure nomenclature scheme ("GROMOS", "AMBER", "CHARMM")
   clean_pdb = True          # Clean the input structure
   pdb2pqr_h_opt = True      # Allow pdb2pqr to optimize the hydrogen positions
   remove_hs = True          # Remove hydrogens in the input structure
   keep_ions = False         # Keep NA+ and CL- ions in the structure
   ser_thr_titration = True  # Titrate SER and THR residues
   
   
   # Poisson-Boltzmann
   gsize = 81           # number of grid nodes
   scaleP = 1           # inverse grid size of big grid
   scaleM = 4           # inverse grid size of small grid
   bndcon = 3           # DelPhi type of boundary condition for the big grid. 3 for Coulombic
   convergence = 0.01   # PB convergence criterion
   nlit = 500           # Number of linear iterations
   nonit = 0            # Number of non-linear iterations
   relfac = 0.75        # DelPhi spectral radius
   relpar = 0.75        # DelPhi relaxation parameter in non-linear iteration convergence process.
   epssol = 80.0        # Solvent dielectric constant value
   
   # Monte Carlo
   pH = 0,14        # pH range: pH minimum, pH maximum OR single pH value
   pHstep = 0.25    # pH range step
   seed = 1234567   # random number generator seed number
   mcsteps = 200000 # number of production monte carlo steps
   eqsteps = 1000   # number of equilibration monte carlo steps
