from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("mc",
              ["mc.pyx"],
              extra_compile_args = ["-fopenmp"],
              extra_link_args=['-fopenmp'],
              include_dirs=[numpy.get_include()]
    )
]

setup(ext_modules=cythonize(extensions))
