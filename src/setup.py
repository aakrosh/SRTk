from distutils.core import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("sr_sam", ["sr_sam.pyx"]),
    Extension("find_best_split", ["find_best_split.pyx"]),
    Extension("combine_alignments", ["combine_alignments.pyx"])
]

setup(version='1.0',  \
      ext_modules = cythonize(extensions),
     )
