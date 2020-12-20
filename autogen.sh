#!/bin/sh
#
# Script for generating the autotools configure script
#

libtoolize --ltdl --force --copy --automake  # version >= 1.5.2
aclocal -I Config                     # version >= 1.11 for Python 3, 1.8.2 for Python 2.6-2.7
autoheader                            # version >= 2.59
automake -a -c                        # version >= 1.11 for Python 3, 1.8.2 for Python 2.6-2.7
autoconf                              # version >= 2.59

# generate bzr revision header file
if bzr version-info
then 
    bzr version-info --custom --template="static int s_bzr_revision={revno};\n" > bzrversion.h
    echo "bzr version header generated" 
elif [ -e bzrversion.h ]
then echo "using existing bzrversion.h"
else 
    echo "static int s_bzr_revision=-1;" > bzrversion.h
    echo "writing dummy bzr revision -1 to header file"
fi
