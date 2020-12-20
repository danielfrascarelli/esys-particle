#!/bin/sh
# A script to move all the output files of a test script into a new folder

mkdir -p $1
mv $2* $1 2>/dev/null
mv snap_* $1 2>/dev/null
