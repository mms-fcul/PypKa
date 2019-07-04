Input Parameters
=================================

============================
General Parameters
============================

.. object:: structure (str)

   :optional: No

   The input filename in Protein Data Bank or GROMACS formats

.. object:: ncpus (int)

   :optional: No

   Number of CPUs to use in the calculations

.. object:: ionicstr (float)

   :optional: No

   Ionic strength of the medium

.. object:: ffinput='GROMOS' (str)

   :allowed values: ('GROMOS', 'AMBER', 'CHARMM')

   Forcefield of the input :py:data:`structure`. GROMOS is
   also valid for experimentally obtained strucutures
		    
.. object:: ffID='G54a7' (str)

   :allowed values: ('G54A7')
      
   For the time being, only GROMOS54a7 based charged and radii are
   allowed.  In the future more forcefields will be added.
		    

.. object:: sites='all' (str)

   CLI Syntax:
   sites_A  = all
   sites_A  = 1N, 18, 35, 48, 66, 129C

   API Syntax:
   sites = 'all'
   sites = {'A': ('1N', '18', '35', '48', '66', '129C')}
   Titration(parameters, sites=sites)

.. object:: clean_pdb=True (boolean)

   The input :py:data:`structure` files are preprocessed with
   PDB2PQR. To skip this step set :py:data:`clean_pdb` to False.
	     
   CLI Syntax: yes or no
   
   
============================
Poisson-Boltzmann Parameters
============================

.. object:: epsin (float)

   :optional: No

   Internal dielectric to be used in the PB calculations.

.. object:: pbc_dimensions (int)

   :optional: No
   :allowed values: (0, 2)

   Number of dimensions with periodic boundaries. Use 0 for solvated
   proteins and 2 for lipidic systems.

.. object:: gsize=81 (int)

   DelPhi number of grid points per side of the cubic lattice.
    
.. object:: epssol=80.0 (float)

   Solvent dielectric to be used in the PB calculations.

.. object:: relpar=0.75 (float)

   DelPhi relaxation parameter in non-linear iteration convergence process.
   
.. object:: relfac=0.75 (float)

   DelPhi spectral radius.
   
.. object:: nonit=0 (float)

   DelPhi number of non-linear iterations.

.. object:: nlit=500 (float)

   DelPhi number of linear iterations.
   
.. object:: convergence=0.01 (float)

   DelPhi convergence threshold for maximum change of potential.
    
.. object:: scaleM=4 (float)

   DelPhi reciprocal of one grid spacing for the finer grid.
    
.. object:: scaleP=1 (float)

   DelPhi reciprocal of one grid spacing for the coarse grid used in the focusing step.
    
.. object:: slice=0.05 (float)

   Only used when :py:data:`pbc_dimensions` is 2. Fraction of the
   exterior of the lipid bilayer to be replicated with periodic
   boundary conditions.
   
.. object:: pbx=False (boolean)

   DelPhi periodic boundary conditions for the x edge of the grid

.. object:: pby=False (boolean)

   DelPhi periodic boundary conditions for the y edge of the grid

.. object:: bndcon=4 (int)

   :allowed values: (1, 2, 3, 4)

   DelPhi type of boundary condition.

   1. Zero
   2. Dipolar
   3. Focusing
   4. Coulombic
		    
.. object:: precision='single' (str)

   :allowed values: ('single', 'double')

   Precision of the compiled DelPhi version being used.


============================
Monte Carlo Parameters
============================
.. object:: seed (int)

   Seed for the MC random generator

.. object:: pH='0,14' (str)

   pH value of the calculation. It can be a single value or a range.
	  
   CLI syntax:
	  
   pH = 0, 14

   
   API syntax:
   
   'pH': '0, 14'


.. object:: pHstep=0.25 (float)

   In case a pH range was provided :py:data:`pH`, the pH step of said
   range is :py:data:`pHstep`
    
