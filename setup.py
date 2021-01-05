from setuptools import setup
from distutils.extension import Extension

setup(
    use_scm_version={
        "write_to": "pypka/_version.py",
        "write_to_template": '__version__ = "{version}"',
        "version_scheme": "post-release",
    },
    packages=["pypka"],
    include_package_data=True,
    setup_requires=["setuptools_scm"],
)