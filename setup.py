#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : setup.py
@Author : XinWang
"""
from setuptools import dist
dist.Distribution().fetch_build_eggs(['Cython', 'numpy'])

import numpy as np
# from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
ext_modules = [Extension("fast_dist",
                         ["fast_dist.pyx"]
                         ),
               Extension('fast_score',
                         ["fast_score.pyx"]
                         )
               ]
setuptools.setup(
    name="ChromoPhyloGen",
    version="1.0",
    author="XinWang",
    author_email="wangx768@mail2.sysu.edu.cn",
    description="A package for inferring CNA fitness evolutionary trees and CNAs' evolutionary efficiency.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FangWang-SYSU/ChromoPhyloGen.git",
    packages=setuptools.find_packages(),
    setup_requires=['Cython', 'numpy'],
    install_requires=['biopython',
                      'Cython',
                      'numpy',
                      'pandas',
                      'scikit-learn',
                      'scipy',
                      'tqdm',
                      'snfpy',
                      'matplotlib',
                      'seaborn',
                      'statsmodels'],
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    include_dirs=[np.get_include()],
    python_requires='>=3.9',
    keywords=['Single cell', 'phylogenetic tree', 'CNA', 'evolutionary efficiency'],
    license='MIT',
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "ChromoPhyloGen=ChromoPhyloGen.cli_main:main",
        ],
    },
)

