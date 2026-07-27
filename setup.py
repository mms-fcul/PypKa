"""Build the Cython Monte Carlo extension.

Everything else is declared in pyproject.toml. The extension is compiled at
build time so that PypKa works on every supported CPython version instead of
relying on pre-built shared objects.
"""

from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
except ImportError:  # pragma: no cover - Cython is a build requirement
    cythonize = None

if cythonize is not None:
    ext_modules = cythonize(
        [Extension("pypka.mc.mc", ["pypka/mc/mc.pyx"])],
        compiler_directives={"language_level": "3"},
    )
else:
    ext_modules = [Extension("pypka.mc.mc", ["pypka/mc/mc.c"])]

setup(ext_modules=ext_modules)
