Main Features
=================================

=================================
Easy to Install
=================================

As seen in the `installation page <installation.html>`_, PypKa is
easilly installed from the command line.

=================================
Easy to Use
=================================

The `basic example <example.html>`_ provided can be effortlessly
adjusted to your system. Although many `parameters <parameters.html>`_
can be modified by an advanced user, simple users can confidently use
the default values.

=================================
API & CLI
=================================

Every interface has its pros and cons. With PypKa you can choose to
run you calculations via a python API or from the command-line. In the
future we will present a webserver interface as well. 

=================================
Fast and accurate
=================================

As a Poisson-Boltzmann-based pKa predictor, PypKa provides an
interesting trade-off between speed and accuracy. You can also
manually set some input parameters such as :py:data:`convergence`,
:py:data:`gsize` or :py:data:`sites`, to adjust the balance between
accuracy and speed as you please. Further speed gains can be achieved
by taking advantage of the highly scalable multiprocessing
capabilities (:py:data:`ncpus`).

=================================
Flexibility of Input Structures
=================================

PypKa supports both the Protein Data Bank and GROMACS input format and
the most popular atomistic force-fields (AMBER, CHARMM & GROMOS) as
well as experimentally determined structures. These structures are
preprocessed using a modified version of PDB2PQR according to the
defined :py:data:`ffinput`.

=================================
Lipidic Systems Support
=================================

It is possible to run calculations on membrane proteins, and in the
future, on lipids.  Currently only some lipids are supported (DMPC,
POPC, POPE and cholesterol) but more will be added. To use this
feature the user is required to name the lipids in their structures
according to the PypKa definition.

Example POPC file

Example POPE file

Example DMPC file

`Example cholesterol file <_static/chol.pdb>`_

