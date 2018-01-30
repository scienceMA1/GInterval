#!/usr/bin/env python

from distutils.core import setup

setup(name='ginterval',
      version='1.0.2',
      description='Genome Interval',
      author='Zonggui Chen',
      author_email='chenzonggui@whu.edu.cn',
      url='https://genome.whu.edu.cn',
      package_dir={'': 'src'},
      py_modules=["ginterval"],
      requires=['pysam', 'biopython']
      )
