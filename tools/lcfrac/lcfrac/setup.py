#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import msnpy


def main():
    
    setuptools.setup(name="lcfrac",
        version=msnpy.__version__,
        description="Python package for processing LC-MS/MS fractionation results",
        long_description=open('README.rst').read(),
        author="Thomas Lawson",
        author_email="t.n.lawson@bham.ac.uk",
        url="https://github.com/computational-metabolomics/lcfrac-galaxy",
        license="GPLv3",
        platforms=['UNIX'],
        keywords=['Metabolomics', 'Lipidomics', 'Mass spectrometry', 'Data Processing',
                  'Direct-Infusion Mass Spectrometry', 'Liquid Chromatography mass spectrometry'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        python_requires='>=3.7',
        install_requires=open('requirements.txt').read().splitlines(),
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
        ],
        entry_points={
         'console_scripts': [
             'lcfrac = lcfrac.__main__:main'
         ]
        }
    )


if __name__ == "__main__":
    main()
