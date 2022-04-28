# ALM

[![License][license-image]][license-url]
[![Doc status][docs-image]][docs-url]
[![Conda build](https://github.com/ttadano/ALM/actions/workflows/build-and-test-ALM-conda.yml/badge.svg)](https://github.com/ttadano/ALM/actions/workflows/build-and-test-ALM-conda.yml)

### Version 2.0.0 Dev

- - -

## Introduction 

This is a software for calculating harmonic and anharmonic interatomic force constants in solids and molecules.

## Features

* Extraction of harmonic and anharmonic force constants based on the supercell approach
* Compressive sensing methods for an efficient and accurate estimation of force constants (LASSO, Elastic net)
* Applicable to any crystal structures and low-dimensional systems
* Accurate treatment of translational and rotational invariance
* Interface to VASP, Quantum-ESPRESSO, OpenMX, xTAPP, and LAMMPS codes
* API for python and C++

## License
Copyright (c) 2014 Terumasa Tadano
This software is released under the MIT license. 
For license rights and limitations, see LICENSE.txt file.

## Author
* Terumasa Tadano
* Atsushi Togo


[license-image]: https://img.shields.io/github/license/ttadano/ALM.svg
[license-url]: https://github.com/ttadano/ALM/blob/develop/LICENSE.txt

[docs-image]: https://readthedocs.org/projects/alm/badge/?version=develop
[docs-url]: https://alm.readthedocs.io/en/develop/?badge=develop

[travis-image]: https://travis-ci.org/ttadano/ALM.svg?branch=develop
[travis-url]: https://travis-ci.org/ttadano/ALM
