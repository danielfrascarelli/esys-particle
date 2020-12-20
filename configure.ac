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

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


# Boilerplate.
AC_INIT([ESyS-Particle], [2.3.5], [d.weatherley@uq.edu.au])
AC_PREREQ([2.59])
AC_LANG([C++])
AC_SUBST([PACKAGE_HOMEPAGE], [https://launchpad.net/esys-particle/])
AC_CONFIG_SRCDIR([tml/type/gettype.cpp])
AC_CONFIG_HEADERS([config.h])

# autoconf generates (in configure) a prefix variable test inside the test for
# the python script directory based on information in python.m4.  Depending on
# the version of this macro definition, the prefix variable test does not pass
# the default prefix value to the prefix variable.  (This happens in the
# python.m4 associated with aclocal 1.9 and 1.10.  The version with aclocal 1.11
# corrects the problem.)  In this case general libraries are installed in
# /usr/lib rather than the default /usr/local/lib, and Python 2.6 libraries are
# installed in /usr/lib/python2.6/dist-packages rather than the default
# /usr/local/lib/python2.6/dist-packages (.../python2.5/site-packages if Python
# 2.5 is being used).
#
# Although autoconf automatically generates a test for the value
# of the prefix variable, this test occurs too late in the configure script if
# based on a problematic python.m4.  Repeating the test here generates the test
# earlier in configure and is a safeguard in the event that a newer python.m4
# returns to the problem encountered with the aclocal 1.10 version.
test "x$prefix" = "xNONE" && prefix=$ac_default_prefix
test "x$exec_prefix" = "xNONE" && exec_prefix=$prefix

# Versions of Automake from at least 1.11 to 1.15 have a bug when using the subdir-objects
# option (see http://debbugs.gnu.org/cgi/bugreport.cgi?bug=13928): A variable appearing in 
# a Makefile.am _SOURCES entry (see Model/Makefile.am for examples) remains unexpanded when 
# running 'configure', being treated instead as a literal directory name and causing the 
# build to fail.  The fix for this bug will appear in Automake 1.16.  When the time arrives
# that Automake requires the subdir-objects option for forward compatibility (perhaps
# version 2.0), hopefully 1.16 or higher will be the default version on operating systems.  
# Failing that, the patch deleted from autogen.sh in revision 1170 can be reinstated to
# remove the subdir-objects option in configure.ac when using Automake 1.15 or older
# (with a small modification for Automake version number, and of course a change to the
# associated comments reflecting the real reason for the fix --- the patch is for build
# error, not for Automake warning).  For now, do not use the option.
#AM_INIT_AUTOMAKE([subdir-objects])
AM_INIT_AUTOMAKE([])
#AM_PATH_CPPUNIT([1.10.0], ,[AC_MSG_ERROR([Dud CppUnit version found, must be >= 1.10.0.])])
AM_PATH_PYTHON([2.6])

#
# These are the new Libtool macros.  Add them and remove the deprecated macros 
# (not forgetting AC_CONFIG_SUBDIRS near the end of this file) if updating
# the Libtool requirements for ESyS-Particle.  Be sure to give the user base
# plenty of warning before requiring a new version of Libtool (which will be
# necessary if the Libtool developers decide no longer to support the deprecated
# macros).
#
#LT_CONFIG_LTDL_DIR([libltdl])
#LT_INIT([dlopen])
#LT_PREREQ([2.0]) 
#LT_LANG([C++])
#LTDL_INIT([subproject convenience])

#
# These Libtool macros are deprecated.  Remove them if updating the Libtool
# requirements for ESyS-Particle.
#
AC_LIBTOOL_DLOPEN
AC_LIBLTDL_CONVENIENCE([libltdl])
AC_PROG_LIBTOOL
AC_SUBST([INCLTDL])
AC_SUBST([LIBLTDL])


#
# Process the escript option --with-escript=<EsysRoot>
#
AC_ARG_WITH(
  [escript],
  AS_HELP_STRING([--with-escript=<EsysRoot>], [Compile with Escript support @<:@default=no@:>@.]),
  [
    case "${withval}" in
      yes) AC_MSG_ERROR([bad value ${withval} for --with-escript, require directory path]) ;;
      no)  withEscript=false; escriptRoot= ;;
      "")  withEscript=true; escriptRoot=. ;;
      *) withEscript=true;escriptRoot=${withval};;
    esac
  ],
  [
    withEscript=false
    escriptRoot=
  ]
)
if test "x$withEscript" = xtrue; then
  #
  # Check whether the path exists.
  #
  AC_CHECK_FILE(
    [$escriptRoot],
    [AC_SUBST(ESCRIPT_ROOT, $escriptRoot)],
    [AC_MSG_ERROR([Path $escriptRoot does not exist!])]
  )
fi
AM_CONDITIONAL([ESCRIPT_ENABLED], [test x$withEscript = xtrue])

#
# Process the VTK option --enable-vtk
#
AC_ARG_ENABLE(
  [vtk],
  AS_HELP_STRING(
    [--enable-vtk],
    [Enable Python Visualization Toolkit support (requires vtk python modules).]
  ),
  [
    if test "x${enableval}" = "xyes"; then
        wants_vtk=true
    else
        wants_vtk=false
    fi
  ],
  [ wants_vtk=false ]
)

if test "x$wants_vtk" = "xtrue" ; then
    AC_PATH_PROG([VTKPYTHON], [vtkpython], "no")
    if test "$VTKPYTHON" = "no" ; then
        AC_MSG_WARN([VTK support enabled but 'vtkpython' executable not found in the path!])
    else
        AC_MSG_NOTICE([VTK support enabled])
    fi
else
    AC_MSG_NOTICE([VTK support disabled])
fi
AM_CONDITIONAL([VTK_ENABLED], [test "x${wants_vtk}" = "xtrue"])

#
# Process the POV-Ray option --disable-povray
#
AC_ARG_ENABLE(
  [povray],
  AS_HELP_STRING(
    [--disable-povray],
    [Disable POV-Ray ray-tracer support @<:@default=enabled if povray found in the PATH@:>@.]
  ),
  [
    if test "x${enableval}" = "xyes"; then
        wants_povray=true
    else
        wants_povray=false
    fi
  ],
  [ wants_povray=true; pov_auto=true ]
)

if test "x$wants_povray" = "xtrue" ; then
    AC_PATH_PROG([POVRAY_EXE], [povray], "no")
    if test "$POVRAY_EXE" = "no" ; then
        if test "x${pov_auto}" = "x" ; then
            AC_MSG_WARN([POV-Ray support enabled but 'povray' executable not found in the path!])
        else
            wants_povray=false
            AC_MSG_NOTICE([POV-Ray support disabled])
        fi
    else
        AC_MSG_NOTICE([POV-Ray support enabled])
    fi
else
    AC_MSG_NOTICE([POV-Ray support disabled])
fi
AM_CONDITIONAL([POVRAY_ENABLED], [test "x${wants_povray}" = "xtrue"])

#
# Process the docs option --enable-docs
#
wants_epydoc=false
wants_pdfdoc=false
AC_ARG_ENABLE(
  [docs],
  AS_HELP_STRING(
    [--enable-docs],
    [Generate HTML and PDF documentation on install (requires epydoc,pdflatex,bibtex)  @<:@default=no@:>@.]
  ),
  [
    if test "x${enableval}" = "xyes"; then
        wants_epydoc=true
        wants_pdfdoc=true
    fi
  ]
)

#
# Process the Epydoc option --with-epydoc[=path/to/epydoc]
#
AC_ARG_WITH(
  [epydoc],
  AS_HELP_STRING(
    [--with-epydoc@<:@=path/to/epydoc@:>@],
    [Generate python API documentation in HTML format on install using epydoc @<:@default=no@:>@.]
  ),
  [
    case "${withval}" in
      no)  wants_epydoc=false ;;
      yes) wants_epydoc=true; epydocDir="" ;;
      *)   wants_epydoc=true; epydocDir="${withval}" ;;
    esac
  ]
)

if test "x$wants_epydoc" = "xtrue" ; then
    if test -n "${epydocDir}" ; then
        AC_PATH_PROG([EPYDOC_EXE], [epydoc], "no", "${epydocDir}")
    else
        AC_PATH_PROG([EPYDOC_EXE], [epydoc], "no")
    fi

    if test "$EPYDOC_EXE" = "no" ; then
        AC_MSG_WARN([epydoc not found! HTML documentation disabled.])
        wants_epydoc=false
    else
        AC_MSG_NOTICE([epydoc found. HTML documentation enabled.])
    fi
else
    AC_MSG_NOTICE([HTML documentation disabled.])
fi
AM_CONDITIONAL([EPYDOC_ENABLED], [test "x${wants_epydoc}" = "xtrue"])

#
# Process the latex option --with-latex[=path/to/latex]
#
AC_ARG_WITH(
  [latex],
  AS_HELP_STRING(
    [--with-latex@<:@=path/to/latex@:>@],
    [Generate documentation in PDF format (requires pdflatex and bibtex) @<:@default=no@:>@.]
  ),
  [
    case "${withval}" in
      yes) wants_pdfdoc=true; latexDir="" ;;
      *)   wants_pdfdoc=true; latexDir="${withval}" ;;
    esac
  ]
)

if test "x$wants_pdfdoc" = "xtrue" ; then
    if test -n "${latexDir}" ; then
        AC_PATH_PROG([PDFLATEX], [pdflatex], "no", "${latexDir}")
        AC_PATH_PROG([BIBTEX], [bibtex], "no", "${latexDir}")
    else
        AC_PATH_PROG([PDFLATEX], [pdflatex], "no")
        AC_PATH_PROG([BIBTEX], [bibtex], "no")
    fi

    if test "$PDFLATEX" = "no" -o "$BIBTEX" = "no"; then
        AC_MSG_WARN([pdflatex and/or bibtex not found! PDF documentation disabled.])
        wants_pdfdoc=false
    else
        AC_MSG_NOTICE([PDF documentation enabled.])
    fi
else
    AC_MSG_NOTICE([PDF documentation disabled.])
fi
AM_CONDITIONAL([PDFDOC_ENABLED], [test "x${wants_pdfdoc}" = "xtrue"])


#
# Process the TAU profiling option --enable-tau=<TauMakefile>
#
AC_ARG_ENABLE(
  [tau],
  AS_HELP_STRING(
    [--enable-tau=<TauMakefile>],
    [Turn on TAU profiling, using TauMakefile to define compiler, compiler-settings, etc.]
  ),
  [
    case "${enableval}" in
      yes) AC_MSG_ERROR([bad value ${enableval} for --enable-tau, require TAU Makefile path]) ;;
      no)  tauProfiling=false; tauMakefile= ;;
      *) tauProfiling=true;tauMakefile=${enableval};;
    esac
  ],
  [
    tauProfiling=false
    tauMakefile=
  ]
)
if test "x$tauProfiling" = xtrue; then
  #
  # Check whether the TAU makefile exists.
  #
  AC_CHECK_FILE(
    [$tauMakefile],
    [AC_SUBST(TAU_MAKEFILE, $tauMakefile)],
    [AC_MSG_ERROR([Makefile $tauMakefile does not exist!])]
  )
fi
AM_CONDITIONAL([TAU_PROFILING], [test x$tauProfiling = xtrue])

#
# Process the SILO option --with-silo
#
silo_prefix=""
AC_ARG_WITH(
	[silo],
  AS_HELP_STRING(
    [--without-silo],
    [Disable SILO database output support @<:@default=detect@:>@]
  ),
  [
    case "${withval}" in
      no) wants_silo=false ;;
      yes) wants_silo=true ;;
      *) wants_silo=true; silo_prefix="${withval}" ;;
    esac
  ],
  [ wants_silo=true ]
)


# Checks for programs.
AC_PROG_CXX([icpc ecpc g++ c++ gpp aCC CC cxx cc++ cl FCC KCC RCC xlC_r xlC])
AC_PROG_CC([icc ecc gcc cc CC cl KCC RCC xlC])
AC_PROG_CPP
AM_PROG_LEX
AC_PROG_YACC
AC_PROG_LN_S
AC_CHECK_TOOL([FIND], [find])
AC_PATH_PROG([LAMINFO], [laminfo])
AM_CONDITIONAL([HAVE_LAM_MPI], [test -n "$LAMINFO"])


#
# Check whether we are using mpic++ or mpiCC as the C++ compiler.
# If we are using mpiCC or mpic++ then don't worry about trying to
# find libmpi and libmpi++
#
AC_MSG_CHECKING(if CXX is mpic++ or mpiCC)
if expr x"$CXX" : x'\(.*mpic[+][+].*\)' '|' x"$CXX" : x'\(.*mpiCC.*\)'; then
  AM_CONDITIONAL([HAVE_SGI_MPI], [false])
else
  AC_MSG_RESULT([...neither: ${CXX}])
  AC_CHECK_LIB([mpi], [MPI_Init], , [AC_MSG_ERROR([MPI library not found...please install MPI])])
  AC_CHECK_LIB([mpi++], [MPI_Abort])
  AC_PATH_PROG([MPIVIS], [mpivis])
  AM_CONDITIONAL([HAVE_SGI_MPI], [test -n "$MPIVIS"])
fi


# Checks for libraries.
PYTHON_SYS_PREFIX=`$PYTHON -c "from __future__ import print_function;import sys;print(sys.prefix)"`
AC_CHECK_LIB([python${PYTHON_VERSION}${PYTHON_VERSION_SUFFIX}], [Py_Main])
PYTHON_CPPFLAGS="-I$PYTHON_SYS_PREFIX/include/python${PYTHON_VERSION}${PYTHON_VERSION_SUFFIX}"

AC_SUBST([PYTHON_CPPFLAGS])

# Check for boost library
if [[ $PYTHON_MAJOR_VERSION -ge 4 ]]; then
  BOOST_REQUIRE([1.47.0])
elif [[ $PYTHON_MAJOR_VERSION -ge 3 ]]; then
  if [[ $PYTHON_MINOR_VERSION -ge 2 ]]; then
    BOOST_REQUIRE([1.47.0])
  else
    BOOST_REQUIRE([1.41.0])
  fi
else
  BOOST_REQUIRE([1.34.1])
fi


# Check for boost filesystem and python libraries
BOOST_FILESYSTEM
BOOST_PYTHON

# check for SILO library
if test "x${wants_silo}" = "xtrue"; then
  CPPFLAGS_SAVED="${CPPFLAGS}"
  LDFLAGS_SAVED="${LDFLAGS}"
  if test "x${silo_prefix}" != "x"; then
    SILO_CPPFLAGS="-I${silo_prefix}/include"
    SILO_LDFLAGS="-L${silo_prefix}/lib"
    AC_SUBST([SILO_CPPFLAGS])
    AC_SUBST([SILO_LDFLAGS])
    CPPFLAGS="${CPPFLAGS} ${SILO_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${SILO_LDFLAGS}"
  fi
  ac_cv_silo=0
  # Try HDF5-enabled SILO first
  AC_CHECK_LIB([siloh5], [DBClose],
    [AC_CHECK_HEADER(silo.h, [ac_cv_silo_h5=1], [ac_cv_silo_h5=0])],
    [ac_cv_silo_h5=0],
    [-lhdf5 -lz]
  )
  if test $ac_cv_silo_h5 -eq 0; then
    # Now check for SILO without HDF5
    AC_CHECK_LIB([silo], [DBClose],
      [AC_CHECK_HEADER(silo.h, [ac_cv_silo=1])]
    )
  fi
  if test $ac_cv_silo_h5 -eq 0 -a $ac_cv_silo -eq 0 ; then
    wants_silo=false
    AC_MSG_WARN([
      The SILO library was not found. SILO output will not be available.
    ])
  else
    if test $ac_cv_silo -eq 1; then
      SILO_LDFLAGS="-lsilo -lm"
    else
      SILO_LDFLAGS="-lsiloh5 -lhdf5 -lz -lm"
    fi
    AC_DEFINE([HAVE_SILO], 1,
      [Define to 1 if the SILO library is available])
  fi
  CPPFLAGS="${CPPFLAGS_SAVED}"
  LDFLAGS="${LDFLAGS_SAVED}"
fi
AM_CONDITIONAL(SILO_ENABLED, test x$wants_silo = xtrue)


# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_HEADER_STDBOOL
AC_HEADER_TIME
AC_CHECK_HEADERS([arpa/inet.h fcntl.h hash_map malloc.h mpi.h netdb.h netinet/in.h stddef.h stdlib.h string.h sys/ioctl.h asm/ioctls.h sys/socket.h sys/time.h unistd.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


# Checks for library functions.
AC_PROG_GCC_TRADITIONAL
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([bzero floor gethostbyname gettimeofday memset pow socket sqrt strdup strerror strtol])
AC_CHECK_FUNCS([MPI_Init], ,[AC_MSG_ERROR([MPI_Init not working...check your MPI installation])])

# We want a versioned documentation directory
lsmdatadocdir="${datarootdir}/doc/${PACKAGE_TARNAME}-${PACKAGE_VERSION}"
AC_SUBST([lsmdatadocdir])

# Output.
AC_CONFIG_SUBDIRS([libltdl]) # Remove when adding the LT* macros above.
AC_CONFIG_FILES([
  Makefile
  ntable/src/Makefile
  ntable/Makefile
  tml/comm/Makefile
  tml/message/Makefile
  tml/type/Makefile
  tml/Makefile
  ppa/src/Makefile
  ppa/Makefile
  pis/Makefile
  Foundation/Makefile
  Geometry/Makefile
  Model/Makefile
  Parallel/Makefile
  Fields/Makefile
  Tools/Makefile
  Tools/dump2geo/Makefile
  Tools/dump2pov/Makefile
  Tools/mesh2pov/Makefile
  Tools/rotextract/Makefile	
  Tools/StressCalculator/Makefile
  Tools/ExtractGrains/Makefile
  Tools/dump2silo/Makefile
  Tools/dump2vtk/Makefile
  Tools/ExtractFractures/Makefile
  Tools/ForceChains/Makefile
  Tools/ExtractStrain/Makefile	
  Tools/MergeLog/Makefile	
  Python/Makefile
  Python/MpiPython/Makefile
  Python/BoostPythonUtil/Makefile
  Python/epydoc/Makefile
  Python/esys/Makefile
  Python/esys/lsm/Makefile
  Python/esys/lsm/util/Makefile
  Python/esys/lsm/geometry/Makefile
  Python/esys/lsm/vis/Makefile
  Python/esys/lsm/vis/core/Makefile
  Python/esys/lsm/vis/vtk/Makefile
  Python/esys/lsm/vis/povray/Makefile
  Doc/Tutorial/Makefile
])
AC_OUTPUT
