#############################################################
##                                                         ##
## Copyright (c) 2003-2017 by The University of Queensland ##
## Centre for Geoscience Computing                         ##
## http://earth.uq.edu.au/centre-geoscience-computing      ##
##                                                         ##
## Primary Business: Brisbane, Queensland, Australia       ##
## Licensed under the Open Software License version 3.0    ##
## http://www.apache.org/licenses/LICENSE-2.0              ##
##                                                         ##
#############################################################
__docformat__ = "restructuredtext en"

import esys.lsm.doc.Util

__installDoc = \
"""
Installation
============

Package Dependencies
--------------------
In order to build the {pkgName:s} software, the following packages
are required:

- boost version 1.32 or later (including boost-python): http://www.boost.org
- python version 2.3 or later: http://www.python.org
- MPI: http://www.mpi-forum.org (known to work with with LAM/MPI version
  6.7 or later http://www.lam-mpi.org and with SGI's MPT
  http://techpubs.sgi.com/)
- CppUnit version 1.10 or later http://cppunit.sourceforge.net/cppunit-wiki

Optional packages include:

- epydoc version 2.1 or later for building python API
  HTML documentation: http://epydoc.sourceforge.net
- docutils version 0.3.5 or later for building tutorial documentation:
  http://docutils.sourceforge.net
- povray version 3.6 or later for rendering simulation data:
  http://www.povray.org
- VTK version 4.2 or later for rendering simulation data:
  http://www.vtk.org

Building the Source Package
---------------------------
{pkgName:s} uses the autotools (autoconf, automake and libtool) build
system to compile source code. After the lsmearth-{version:s}.tar.gz source
has been downloaded and unziped/untared, change into the
{pkgName:s} top level source directory which contains the ``configure``
script. The ``configure`` shell script attempts to guess correct values for
various system-dependent variables used during compilation.
A typical invocation of the ``configure`` script on a system with LAM-MPI
installed would be::

  $ ./configure CC=mpicc CXX=mpic++ --prefix=/opt/pkgs/lsmearth-{version:s}

On an MPT system (Altix 3700, for example), the invocation of ``configure``
would be::

$ ./configure CC=gcc CXX=g++ --prefix=/opt/pkgs/lsmearth-{version:s}

In both these cases, the ``CC`` and ``CXX`` configure variables are used to set the
C and C++ compilers, respectively. The ``--prefix=/opt/pkgs/lsmearth-{version:s}``
option sets the root installation path, so that in these cases the {pkgName:s}
executables, libraries, etc will be installed under ``/opt/pkgs/lsmearth-{version:s}/bin``,
``/opt/pkgs/lsmearth-{version:s}/lib``, etc.
Running::

$ ./configure --help

will list the command line options (and brief descriptions) for the ``configure`` script.

Once the ``configure`` script has been run, and the various ``Makefile`` files
have been generated, the binaries can be built by issuing the::

$ make

command. Parallel builds are supported, so on multiple (or hyper-threaded)
CPU machines::

$ make -j N

will cause ``N`` files to compile simultaneously. Issuing the::

$ make install

command installs binary, data and python package and module files.
The installation can be tested with::

$ make installcheck

which runs a series of unit-tests using the **installed** package binaries.

Environment Variables
---------------------
Typically, there are three environment variables which need to
be set before being able to run an {pkgName:s} python script: ``PATH``,
``LD_LIBRARAY_PATH`` and ``PYTHONPATH``. Assuming that the ``configure``
script was run with the ``--prefix=/opt/pkgs/lsmearth-{version:s}``
option, these environment variables would be set as follows::

$ export PATH=/opt/pkgs/lsmearth-{version:s}/bin:$PATH
$ export LD_LIBRARAY_PATH=/opt/pkgs/lsmearth-{version:s}/lib:$LD_LIBRARY_PATH
$ export PYTHONPATH=/opt/pkgs/lsmearth-{version:s}/lib/python2.3/site-packages:$PYTHONPATH

:summary: Installation documentation.

"""

__doc__ = esys.lsm.doc.Util.setSectionDoc("InstallSection", __installDoc)
