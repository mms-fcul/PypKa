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

* python3.9 up to python3.14
* libgfortran4
* gawk
* a C compiler (only required when building from source)

These can all be installed via apt in Debian based systems::

  apt install gawk gcc gfortran libgfortran4

Via Docker
----------

Run the following command on a local directory which contains the cli parameters file (<PARAMETERS_FILE>) and input structure::

   docker run -v ${PWD}:/home/ -w /home -t pedrishi/pypka:latest python3 -m pypka <PARAMETERS_FILE>
