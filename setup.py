from distutils.core import setup, Extension
import numpy

setup(name = 'pyGravSphere', author = 'Anna Genina, Justin Read', author_email = 'anna.genina@durham.ac.uk', ext_modules=[Extension("_gravsphere",
      sources=["gravsphere.c", "gravsphere.i"],
      include_dirs=[numpy.get_include()])])
