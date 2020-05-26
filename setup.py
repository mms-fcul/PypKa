import setuptools

long_description = """
# PypKa

A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism


# Dependencies

  - libgfortran4
  - gawk
  - pytest
  - numpy
  - psutil


# License

  pypka is distributed under a LGPL-3.0, however delphi4py depends on
  DelPhi which is proprietary. To use DelPhi the user is required to
  download the DelPhi license
  [here](https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi)

# Documentation

  Documentation can be found [here](https://pypka.readthedocs.io/en/latest/). (Under development)

# Installation

  pip3 install pypka

  apt install gawk gcc gfortran libgfortran4

# Contacts

  Please submit a github issue to report bugs and to request new features.
  Alternatively you may find the developer [here](mailto:pdreis@fc.ul.pt). Please visit ou [website](http://mms.rd.ciencias.ulisboa.pt/) for more information.

"""



setuptools.setup(
    name="pypka",
    version="1.2.0",
    author="Pedro Reis",
    author_email="pdreis@fc.ul.pt",
    description="A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mms-fcul/PypKa",
    packages=setuptools.find_packages(),
    package_dir={'pypka': 'pypka'},
    package_data={'pypka': ['mc*.so',
                            'addHtaut*',
                            'G54A7/*',
                            'G54A7/sts/*',
                            'pdb2pqr/*',
                            'pdb2pqr/dat/*',
                            'pdb2pqr/extensions/*',
                            'pdb2pqr/src/*',
                            'pdb2pqr/ZSI/*',
                            'delphi4py/*',
                            'delphi4py/readFiles/*.so',
                            'delphi4py/rundelphi/*.so']},
    install_requires=['numpy', 'pytest', 'coverage', 'psutil'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3) ",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    entry_points={
        'console_scripts': [
            'pypka = pypka.pypka:CLI'
        ]
    }
)
