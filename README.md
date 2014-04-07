ADANET algorithm
======

This is an implementation of ADANET algorithm for Gene Regulatory Network inference from mRNA expression data, in form of a set of Matlab scripts and C++ workers. See: http://dl.acm.org/citation.cfm?id=2382992 for more information.

Requirements:
- Matlab environment, the script was implemented and tested using Matlab 7.10.0.499 (R2010a).
- c++ compiler and standard libraries


Installation:
- Download and unpack the source code to the download folder.
- Build the adanet mex file by invoking the appropriate command in the main adanet folder, i.e. a folder containing Makefile, from command line:

$make

Usage:
- Once the script is in working path, please refer to the example.m file.

LICENSE
=====

Copyright (C) 2014  Janusz Slawek

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program, see LICENSE.
