from distutils.core import setup
from Cython.Build import cythonize

from pore.paths import ROOT_DIR

setup(ext_modules = cythonize(str(ROOT_DIR / "pore" / "cython" / "cython_voxel.pyx")))