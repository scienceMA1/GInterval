#!/usr/bin/env python

import setuptools

setup_params = dict(
    name="ginterval",
    version="1.0.3",
    #description="Genomic interval",
    author="Zonggui Chen",
    #author_email="chenzonggui@whu.edu.cn",
    #keywords="omit",
    #url="omit",
    package_dir={'': 'src'},
    #src_root=None,
    packages=setuptools.find_packages("src"),
    include_package_data=True,
    test_suite="nose.collector"
    #py_packages=["ginterval"],
)


if __name__ == '__main__':
    dist = setuptools.setup(**setup_params)

