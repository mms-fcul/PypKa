import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pypka",
    version="0.0.4",
    author="Pedro Reis",
    author_email="pdreis@fc.ul.pt",
    description="A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mms-fcul/PypKa",
    packages=setuptools.find_packages(),
    package_dir={'pypka': 'pypka'},
    package_data={'pypka': ['mc.so',
                            'addHtaut',
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
    install_requires=['numpy', 'pytest', 'coverage'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
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
