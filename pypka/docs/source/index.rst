.. PypKa documentation master file, created by
   sphinx-quickstart on Tue Mar 26 14:52:09 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Content
   :hidden:

   self
   installation
   example
   features
   methods
   parameters
   future
   GitHub Project <https://github.com/mms-fcul/PypKa>

Pypka
=================================

A python module for flexible Poisson-Boltzmann based p\ :emphasis:`K`\
:sub:`a` calculations.

We have implemented a flexible tool to predict Poisson-Boltzmann-based
p\ :emphasis:`K`\ :sub:`a` values of biomolecules. This is a free and
open source project that provides a simple, reusable and extensible
python API and CLI for p\ :emphasis:`K`\ :sub:`a` calculations with a
valuable trade-off between fast and accurate predictions. With PypKa
one can enable p\ :emphasis:`K`\ :sub:`a` calculations, including
optional proton tautomerism, within existing protocols by adding a few
extra lines of code. PypKa supports CPU parallel computing on
anisotropic (membrane) and isotropic (protein) systems and allows the
user to find a balance between accuracy and speed.

PypKa is written in Python and Cython and it is integrated with the
Poisson-Boltzmann solver DelPhi Fortran77 via `DelPhi4Py
<https://github.com/mms-fcul/DelPhi4Py>`_.


=================
Availability
=================

PypKa can be `easily installed <installation.html>`_ using the pip
package manager.


Source code is freely available at `GitHub <https://github.com/mms-fcul/PypKa>`_
under the LGPL-3.0 license. The package can be installed from `PyPi
<https://pypi.org/project/pypka/>`_.

=================
Contacts
=================

Please `submit a github issue  <https://github.com/mms-fcul/PypKa/issues/new>`_ to report bugs and to request new features. 

Alternatively you may find the developer `here <mailto:pdreis@fc.ul.pt>`_ . Please visit our `website <http://mms.rd.ciencias.ulisboa.pt/>`_ for more information.

