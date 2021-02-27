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
   newresidues
   contributing
   future
   JCIM Article <https://doi.org/10.1021/acs.jcim.0c00718>
   GitHub Project <https://github.com/mms-fcul/PypKa>
   pypka

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

If you use PypKa in your research please cite the following `paper <https://doi.org/10.1021/acs.jcim.0c00718>`_:

Reis, P. B. P. S., Vila-Viçosa, D., Rocchia, W., & Machuqueiro, M. (2020). 
PypKa: A Flexible Python Module for Poisson–Boltzmann-Based pKa Calculations. 
Journal of Chemical Information and Modeling, 60(10), 4442–4448.

.. code-block:: bibtex

   @article{reis2020jcim,
   author = {Reis, Pedro B. P. S. and Vila-Viçosa, Diogo and Rocchia, Walter and Machuqueiro, Miguel},
   title = {PypKa: A Flexible Python Module for Poisson–Boltzmann-Based pKa Calculations},
   journal = {Journal of Chemical Information and Modeling},
   volume = {60},
   number = {10},
   pages = {4442-4448},
   year = {2020},
   doi = {10.1021/acs.jcim.0c00718}
   }



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

