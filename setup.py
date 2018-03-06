#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import dirname, realpath, exists
from setuptools import setup, find_packages
import sys


author = u"Paul Müller"
authors = [author, "Gheorghe Cojoc"]
description = 'global geometric factors and corresponding stresses of the optical stretcher '
name = 'ggf'
year = "2018"

sys.path.insert(0, realpath(dirname(__file__))+"/"+name)
from _version import version

if __name__ == "__main__":
    setup(
        name=name,
        author=author,
        author_email='dev@craban.de',
        url='https://github.com/GuckLab/ggf',
        version=version,
        packages=find_packages(),
        package_dir={name: name},
        include_package_data=True,
        license="MIT",
        description=description,
        long_description=open('README.rst').read() if exists('README.rst') else '',
        install_requires=["numpy >= 1.9.0",
                          "scipy >= 0.18.0",
                          ],
        setup_requires=['pytest-runner'],
        tests_require=["pytest"],
        python_requires='>=3.5, <4',
        keywords=["optical stretcher",
                  "global geometric factor",
                  ],
        classifiers= [
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Intended Audience :: Science/Research'
                     ],
        platforms=['ALL'],
        )
