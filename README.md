[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mms-fcul/PypKa/blob/master/examples/notebook/pypka.ipynb) [![CircleCI](https://circleci.com/gh/mms-fcul/PypKa.svg?style=svg)](https://circleci.com/gh/mms-fcul/PypKa) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/59a058e4bf0846f18d9d1f6b16a4a0e5)](https://www.codacy.com/gh/mms-fcul/PypKa/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mms-fcul/PypKa&amp;utm_campaign=Badge_Grade) [![Documentation Status](https://readthedocs.org/projects/pypka/badge/?version=latest)](https://pypka.readthedocs.io/en/latest/?badge=latest)

[![PyPI version](https://badge.fury.io/py/pypka.svg)](https://badge.fury.io/py/pypka)  [![PyPI - Downloads](https://img.shields.io/pypi/dm/pypKa)](https://badge.fury.io/py/pypKa)

# PypKa

A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism
DOI: <a href="https://doi.org/10.1021/acs.jcim.0c00718">10.1021/acs.jcim.0c00718</a>

```bibtex
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
```

## Documentation & Basic Usage

  Documentation can be found [here](https://pypka.readthedocs.io/en/latest/)

  Starting templates for the the API and CLI usage can be found [here](https://pypka.readthedocs.io/en/latest/example.html) while a online notebook is hosted at [Google Colab](https://colab.research.google.com/github/mms-fcul/PypKa/blob/master/pypka/example/notebook/pypka.ipynb)

## Installation & Dependencies

PypKa should be installed in a Linux-based system. In Windows please use the [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

* python2.6>= & python3.5>=
* libgfortran4
* gawk
* pytest
* numpy

```bash
  python3 -m pip install pypka
  apt install gawk gcc gfortran libgfortran4 python2
```

A docker image is also [available](https://hub.docker.com/r/pedrishi/pypka). 
A simple way of using it would be to build the image from the [Dockerfile](./Dockerfile) and running it on a local directory which contains the cli parameters file and input structure.
```bash
docker build -t pypka/latest -f Dockerfile .
docker run -v ${PWD}:/home/ -w /home -t pypka:latest python3.9 -m pypka <PARAMETERS_FILE>
```

Functioning working examples of the API, CLI, docker and notebook can be found under /examples.

## Contributing

Contributions are encouraged, and they are greatly appreciated!

You can contribute in many ways:

* Report Bugs
* Fix Bugs
* Implement Features
* Write Documentation
* Submit Feedback

For more info check [CONTRIBUTING](./CONTRIBUTING.rst)

## License

  pypka is distributed under a [LGPL-3.0](./LICENSE), however delphi4py depends on
  DelPhi which is proprietary. To use DelPhi the user is required to
  download the DelPhi license
  [here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

## Contacts

  Please submit a github issue to report bugs and to request new features.
  Alternatively you may find the developer [here](mailto:pdreis@fc.ul.pt). Please visit ou [website](http://mms.rd.ciencias.ulisboa.pt/) for more information.
