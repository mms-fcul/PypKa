Installation
=================================

====================
Python Package Index
====================

To install the latest stable release with pip do::

   pip3 install pypka

====================
Required Software
====================

Pypka depends on the following software:

* Python3.6
* Python2
* gawk

For the moment DelPhi depends on the following software:

* gcc
* gfortran
* libgfortran4

These can all be installed via apt in Debian based systems::

  apt install gawk gcc gfortran libgfortran4 python2

====================
Cross Platform
====================

PypKa has been developed for Linux users and it does not support Windows or MacOS. However, we do provide a docker image for these users.

   docker container run -it -v ${PWD}:/pypka pedrishi/pypka:latest
