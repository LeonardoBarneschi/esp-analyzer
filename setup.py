from setuptools import setup, find_packages

setup(name='grid_analyzer',
      version='0.1',
      description='Perform a number of useful grid operations',
      url='http://github.com/storborg/funniest',
      author='Leonardo Barneschi',
      author_email='leonardo.barneschi@student.unisi.it',
      license='MIT',
      packages=find_packages(),
      scripts=['grid_analysis/bin/grid_slicer',
               'grid_analysis/bin/grid_esp',
               'grid_analysis/bin/esp_3d_plot',
               'grid_analysis/bin/vdw_plane',
               'grid_analysis/bin/computephi'],
      zip_safe=False)
