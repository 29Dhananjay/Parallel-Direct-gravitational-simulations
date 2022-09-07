import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "grav",
        ["OpenMP.pyx"],
        extra_compile_args=['-fopenmp','-v','-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include'],
        extra_link_args=['-lgomp', '-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/10/'],
        include_dirs=[numpy.get_include()],
    )
]

setup(name="OpenMP",
      ext_modules=cythonize(ext_modules))

