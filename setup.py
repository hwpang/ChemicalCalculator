#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="ChemicalCalculator",
    version="0.0.1",
    author="Hao-Wei Pang",
    author_email="hwpang@mit.com",
    description="A chemical calculator for organic peroxide",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hwpang/chemical_calculator.git",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
)
