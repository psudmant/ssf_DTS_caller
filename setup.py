from distutils.core import setup, Extension 
from Cython.Distutils import build_ext
from Cython.Build import cythonize 

ext_modules = [
               Extension("traverse_contours", 
                         ["traverse_contours.pyx"],
                         language="c"),
               Extension("get_windowed_variance", 
                         ["get_windowed_variance.pyx"],
                         language="c"),
               Extension("c_hierarchical_edge_merge",
                         ['c_hierarchical_edge_merge.pyx',
                         './heap/heap.cc',
                         './edge.cc',
                         './float_heap.cc'],
                        language = "c++")
              ]

setup( 
    cmdclass = {'build_ext':build_ext},
    name = 'c_SSF',
    ext_modules = ext_modules 
    )
