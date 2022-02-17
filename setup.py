from setuptools import setup, find_packages
from numpy.distutils.core import Extension, setup

setup(name='grid_analyzer',
      version='0.1',
      description='Perform a number of useful grid operations',
      url='http://github.com/storborg/funniest',
      author='Leonardo Barneschi',
      author_email='leonardo.barneschi@student.unisi.it',
      license='MIT',
      packages=find_packages(),
      ext_modules=[ Extension('grid_analysis.geom_operations.getpot', ['grid_analysis/geom_operations/getpot.f90']),
                    Extension('grid_analysis.geom_operations.getpot_parallel', ['grid_analysis/geom_operations/getpot_parallel.f90'],
                    extra_f90_compile_args=['-fopenmp', '-lgomp'],
                    extra_link_args=['-lgomp']),
                    Extension('grid_analysis.geom_operations.getcoulene_parallel', ['grid_analysis/geom_operations/getcoulene_parallel.f90'],
                    extra_f90_compile_args=['-fopenmp', '-lgomp'],
                    extra_link_args=['-lgomp']),
                    Extension('grid_analysis.geom_operations.getcoulene', ['grid_analysis/geom_operations/getcoulene.f90']),
                    Extension('grid_analysis.file_operations.getmap', ['grid_analysis/file_operations/getmap.f90']),
                    ],

      scripts=['grid_analysis/bin/grid_slicer',
               'grid_analysis/bin/grid_esp',
               'grid_analysis/bin/esp_3d_plot',
               'grid_analysis/bin/vdw_plane',
               'grid_analysis/bin/computephi'],
      zip_safe=False)
