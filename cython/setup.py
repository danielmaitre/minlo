from distutils.core import setup, Extension
from Cython.Build import cythonize

ext=Extension("sudakovs",
              sources=["sudakovs.pyx","../sudakovs.cpp",],
              include_dirs=["/home/daniel/Tools/include"],
              libraries=['MINLOlib','LHAPDF','gsl','gslcblas'],
              extra_link_args=['-L/home/daniel/Tools/lib/','-L/home/daniel/workspace/minlo/build/install_dir/lib/','-Wl,-rpath','-Wl,/home/daniel/workspace/minlo/build/install_dir/lib/']
)

setup(
    name="sudakovs_cpp",
    ext_modules=cythonize(
        [ext],

    )
)
              
