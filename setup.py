from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

# get the annotated file as well
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

ext_modules = [
    Extension('_decomphelmholtz',
              sources            = ['_decomphelmholtz.pyx'],
              include_dirs       = [numpy.get_include(),'.'],
              extra_compile_args = ['-std=c99','-fopenmp','-pthread','-fPIC','-mtune=native','-march=native','-O3','-falign-functions=64'],
              extra_link_args    = ['-fopenmp','-pthread'],
              libraries          = ['fftw3','gomp'],
              library_dirs       = ['/usr/local/lib']),
    
]
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

import os
import shutil

srcfile = 'build/lib.linux-x86_64-2.7/_decomphelmholtz.so'
dstfile = './_decomphelmholtz.so'


assert not os.path.isabs(srcfile)
shutil.copy(srcfile, dstfile)
shutil.rmtree('build')
