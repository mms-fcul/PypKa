Installation
============

There are three ways to install pypka:

- with pip. Preferred for Linux & Windows (with WSL) users.

- with Docker.

- with no installation, by using the online server at `pypka.org <https://pypka.org>`_. Recommended when only a few calculations are going to be performed.




Via pip
-------

To install the latest stable release with pip do::

   pip3 install pypka


Required Software
~~~~~~~~~~~~~~~~~

Pypka depends on the following software:

* python2.6>= & python3.5>=
* libgfortran4
* gawk

These can all be installed via apt in Debian based systems::

  apt install gawk gcc gfortran libgfortran4 python2

Via Docker
----------

Run the following command on a local directory which contains the cli parameters file (<PARAMETERS_FILE>) and input structure::

   docker run -v ${PWD}:/home/ -w /home -t pedrishi/pypka:latest python3.9 -m pypka <PARAMETERS_FILE>
