from distutils.extension import Extension

from setuptools import setup

setup(
    packages=["pypka"],
    include_package_data=True,
    setup_requires=["setuptools_scm"],
    version='2.1.1',
)