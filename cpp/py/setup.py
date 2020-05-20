#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#radia = Extension(
#    'radia',
#    define_macros=[('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
#    #include_dirs=[os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse')],
#    #include_dirs=[os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse'), os.path.abspath('../src/core')], #OC03112019 requested by R. Nagler
#    include_dirs=[os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse'), os.path.abspath('../src/core'), os.path.abspath('/usr/lib/openmpi/include'), os.path.abspath('/usr/lib/openmpi/include/openmpi')], #MPI comp. test #OC03112019 requested by R. Nagler
#    #libraries=['radia', 'm', 'fftw'],
#    libraries=['radia', 'm', 'fftw', 'mpi_cxx', 'dl', 'hwloc'], #MPI comp. test #OCTEST 12/012020
#    library_dirs=[os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')],
#    #library_dirs=[os.path.abspath('../gcc'), os.path.abspath('../../ext_lib'), os.path.abspath('/usr/lib/openmpi/lib')], #MPI comp. test
#    sources=[os.path.abspath('../src/clients/python/radpy.cpp')],
#    extra_compile_args=['-v', '--verbose']) #MPI comp. test 

ext_kwargs = {'define_macros': [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
              #'include_dirs': [os.path.abspath('../src/lib')],
              'include_dirs': [os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse'), os.path.abspath('../src/core')], #os.path.abspath('/usr/lib/openmpi/include'), os.path.abspath('/usr/lib/openmpi/include/openmpi')], #MPI comp. test #OC03112019 requested by R. Nagler
              'libraries': ['radia', 'm', 'fftw'],
              'library_dirs': [os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')],
              'sources': [os.path.abspath('../src/clients/python/radpy.cpp')]} 

if 'MODE' in os.environ: 
    sMode = str(os.environ['MODE'])
    if sMode == 'mpi': 
        ext_kwargs.update({
            'include_dirs': [os.path.abspath('../src/lib'), os.path.abspath('../src/ext/auxparse'), os.path.abspath('../src/core'), os.path.abspath('/usr/lib/openmpi/include')], #os.path.abspath('/usr/lib/openmpi/include/openmpi')], #MPI comp. test #OC03112019 requested by R. Nagler
            'libraries': ['radia', 'm', 'fftw', 'mpi_cxx', 'dl'], #for compilation with OpenMPI, tested on servers at BNL
            'library_dirs': [os.path.abspath('../gcc'), os.path.abspath('../../ext_lib'), os.path.abspath('/usr/lib/openmpi/lib')]})
    elif sMode == 'mpi_nersc': #for compilation with MPICH, tested on NERSC cluster
        ext_kwargs.update({
            'libraries': ['radia', 'm', 'fftw', 'mpichcxx_intel', 'dl'],
            'library_dirs': [os.path.abspath('../gcc'), os.path.abspath('../../ext_lib'), os.path.abspath(os.getenv('MPICH_DIR') + '/lib')]})
    elif sMode == '0':
        pass
        #ext_kwargs.update({'libraries': ['srw', 'm', 'fftw3f', 'fftw3']}) #OC07022019
    else:
        raise Exception("Unknown Radia compilation/linking option")

radia = Extension('radia', **ext_kwargs)

setup(name='Radia Python Interface',
      version='1.0',
      description='This is Radia for Python',
      author='O. Chubar, P. Elleaume, J. Chavanne',
      author_email='chubar@bnl.gov',
      url='http://github.com/ochubar/Radia',
      long_description='''
This is Python interface to the Radia 3D magnetostatic code.
''',
      ext_modules=[radia])
