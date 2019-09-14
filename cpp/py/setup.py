from distutils.core import setup, Extension
import os

radia = Extension(
    'radia',
    define_macros=[('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
    include_dirs=[os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse')],
    libraries=['radia', 'm', 'fftw'],
    library_dirs=[os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')],
    sources=[os.path.abspath('../src/clients/python/radpy.cpp')])

setup(name='Radia Python Interface',
      version='1.0',
      description='This is Radia for Python',
      author='O. Chubar, P. Elleaume, J. Chavanne',
      author_email='chubar@bnl.gov',
      url='http://github.com/ochubar/Radia',
      long_description='''
This is Radia for Python.
''',
      ext_modules=[radia])
