from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = "mc",
  cmdclass = {"build_ext": build_ext},
  ext_modules =
  [
    Extension("mc",
              ["mc.pyx"]
              #extra_compile_args = ["-fopenmp"],
              #extra_link_args=['-fopenmp']
              )
  ]
)
