#!/bin/sh
# generate bzr revision header file
cd $1

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
