from setuptools import setup

setup(
    use_scm_version={
        "write_to": "version.py",
        "write_to_template": '__version__ = "{version}"',
        "version_scheme": "post-release",
    },
    setup_requires=["setuptools_scm"],
)