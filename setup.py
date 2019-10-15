from distutils.core import setup, Extension
import numpy

setup(ext_modules=[Extension("_gravsphere",
      sources=["gravsphere.c", "gravsphere.i"],
      include_dirs=[numpy.get_include()])])
